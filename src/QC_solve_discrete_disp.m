function [disp_n1,result,QCbrokenLink]=QC_solve_discrete_disp(QC_projectCfg,QCboundaryCfg,computeCfg,materialCfg,...
    QCelemCfg,QCnodeCoord,QCelemNodeNo,QCelemNodePos,QCelemDofLabel,QCelemLength,disp_n,N,pNpxi,pNpeta,pNpzeta,M,assem,display,QCbrokenLink)
if display
    disp('Solving QC displacement field.');
end
tmpBrokenLink=QCbrokenLink;
lLc=zeros(QCelemCfg.nElem,1);
brokenElemExpandLimit=[1,1,2,3,4,5,5,6];
Voigt=[1,1;2,2;1,2];
meshDispDOF=QCelemCfg.meshDispDOF;
elemDispDOF=QCelemCfg.elemDispDOF;
DOF=QCelemCfg.phyDOF;
nElem=QCelemCfg.nElem;
nNode=QCelemCfg.nNode;
nElemNode=QCelemCfg.nElemNode;
nGaussPoint=QCelemCfg.nGaussPoint;
W=QCelemCfg.gaussWeight;

kB=materialCfg.kB;
T=materialCfg.T;
b=materialCfg.b;
Lc=materialCfg.Lc;
Lc_L=materialCfg.Lc_L;

% set displacement boundary
[disp_n,saveConsValue]=set_disp_boundary(disp_n,QCboundaryCfg,computeCfg);

disp_init=disp_n;
newtoniter = 0;
while(newtoniter<=10)
    % Initiation
    %     deformed_nodeCoordRef=nodeCoordRef+disp_n;
    %     elemNodePos=deformed_nodeCoordRef(elemDofLabel);
    newtoniter=newtoniter+1;
    Kdisp_quad=zeros(elemDispDOF^2,nElem);
    Rdisp_quad=zeros(elemDispDOF,nElem);
    %%%%%%%%%%%%% optional %%%%%%%%%%%%%%
    fnet_quad=zeros(length(QCelemNodeNo),1);
    H_quad=zeros(nGaussPoint,nElem);
    Hneg_quad=zeros(nGaussPoint,nElem);
    %%%%%%%%%%% end optional %%%%%%%%%%%%
    disp_net=[disp_n(1:2:end),disp_n(2:2:end)];

    parfor iElem=1:nElem
        if ~(QC_projectCfg.fracture && QCbrokenLink(iElem)==1)
            material=cell2mat(materialCfg.materials(1,1));
            Kdisp_e=zeros(elemDispDOF);
            Rdisp_e=zeros(elemDispDOF,1);
            %%%%%%%%%%%%%%%%%%%%%% QC %%%%%%%%%%%%%%%%%%%%%%%%
            netElemNodeNo=QCelemNodeNo(iElem,:);
            netElemDisp=[disp_net(netElemNodeNo(:),:)];
            netElemNodePos=QCnodeCoord(netElemNodeNo,:);
            deformed=netElemNodePos+netElemDisp;
            vec=deformed(1,:)-deformed(2,:);
            origin_vec=netElemNodePos(1,:)-netElemNodePos(2,:);
            l=norm(vec);
            L=QCelemLength(iElem);
            Ncurr=vec/norm(vec);
            Norig=origin_vec/norm(origin_vec);
            lbd=l/L;
            x=l-L;
            Lc=L*Lc_L;
            if abs(x/Lc-1)<0.01
                x=0.9*Lc;
            end
            f=kB*T/b*(1/4*(1-x/Lc)^(-2)-1/4+x/Lc);
            fnet_quad(iElem,1)=f*b/kB/T;
            k=kB*T/b*(1/2/Lc*(1-x/Lc)^(-3)+1/Lc);

            if newtoniter>1
                lLc(iElem)=max(lLc(iElem),l/Lc);
            end
            S=zeros(2,2);
            Cvoigt=zeros(3,3);
            if (QCbrokenLink(iElem)==0 && l/Lc < QC_projectCfg.critlbd) || ~QC_projectCfg.fracture

                S=S+1/lbd*f*L*kron(Norig',Norig);
                for p=1:3
                    for q=1:3
                        if p>=q
                            i=Voigt(p,1);
                            j=Voigt(p,2);
                            m=Voigt(q,1);
                            n=Voigt(q,2);
                            Cvoigt(p,q)=Cvoigt(p,q)+(1/lbd^2*k*L^2-1/lbd^3*f*L)*Norig(i)*Norig(j)*Norig(m)*Norig(n);
                        end
                    end
                end
                for p=1:3
                    for q=1:3
                        if p<q
                            Cvoigt(p,q)=Cvoigt(q,p);
                        end
                    end
                end
                Svoigt=[S(1,1);S(2,2);S(1,2)];

                Bcal=zeros(2,2);
                B0=zeros(3,4);
                x1=netElemNodePos(1,1); x2=netElemNodePos(2,1);
                y1=netElemNodePos(1,2); y2=netElemNodePos(2,2);

                %         N1=@(x,y)1-sqrt((x-x1)^2+(y-y1)^2)/L;
                %         N2=@(x,y)sqrt((x-x1)^2+(y-y1)^2)/L;
                %         pN1px=@(x,y)(x-x1)/sqrt((x-x1)^2+(y-y1)^2)/L;
                %         pN1py=@(x,y)(y-y1)/sqrt((x-x1)^2+(y-y1)^2)/L;

                %         sinT=(y2-y1)/L;
                %         cosT=(x2-x1)/L;
                %         F=lbd*kron(Ncurr',Norig)+0*[sinT^2,-sinT*cosT;-sinT*cosT,cosT^2];
                F=lbd*kron(Ncurr',Norig);
                lx=(x2-x1)/2; ly=(y2-y1)/2;
                pN1px=-lx*2/L^2;
                pN1py=-ly*2/L^2;
                pN2px=-pN1px;
                pN2py=-pN1py;
                pNpX=[pN1px,pN2px;pN1py,pN2py];
                Bu=zeros(3,QCelemCfg.elemDispDOF);
                B0=zeros(3,QCelemCfg.elemDispDOF);
                for iNode=1:QCelemCfg.nElemNode
                    tmp=[pNpX(1,iNode),0;...
                        0,pNpX(2,iNode);...
                        pNpX(2,iNode),pNpX(1,iNode)];
                    Bu(:,iNode*2-1:iNode*2)=tmp;
                    B0(:,iNode*2-1:iNode*2)=tmp*F';
                end

                Bcal=pNpX;

                detJ=1;

                iGauss=1;
                Kgeo=zeros(DOF*nElemNode,DOF*nElemNode);
                for iC=1:nElemNode
                    for jC=1:nElemNode
                        Kij=eye(DOF,DOF);
                        Kij=Kij*(W(iGauss)*Bcal(:,iC)'*S*Bcal(:,jC)*detJ);
                        Kgeo((DOF*iC-(DOF-1)):(DOF*iC),(DOF*jC-(DOF-1)):(DOF*jC))=Kij;
                    end
                end
                Kdisp_e=Kdisp_e+W(iGauss)*B0'*Cvoigt*B0*detJ+Kgeo;
                Rdisp_e=Rdisp_e+W(iGauss)*B0'*Svoigt*detJ;

                Kdisp_quad(:,iElem)=Kdisp_e(:);
                Rdisp_quad(:,iElem)=Rdisp_e;
            elseif newtoniter>1
                tmpBrokenLink(iElem)=1;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%% End QC %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% Assembling %%%%%%%%%%%%%%%%%%%%%%%%%
    Kdisp=sparse(assem.m8_8i(:),assem.m8_8j(:),Kdisp_quad(:),meshDispDOF,meshDispDOF);
    Rdisp=sparse(assem.m8_1(:),1,Rdisp_quad(:),meshDispDOF,1);
    % Interpolating stress from Gauss points to nodes, the denominator depends
    % on the number of Gauss points for interpolation of each point.
    result.fnet=fnet_quad;
    result.reaction=full(Rdisp);
    %%%%%%%%%%%%%%%%% end assembling %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Newton-Raphson %%%%%%%%
    brokenElems=find(QCbrokenLink);
    delDOF=QCboundaryCfg.dispConstraint;
    brokenDOF=[];
    saveBrokenValue=[];
    if QC_projectCfg.fracture && ~isempty(brokenElems)
        forceBroken=false;
        influencedNodes=QCelemNodeNo(brokenElems,:);
        nodesCount=zeros(QCelemCfg.nNode,1);
        tab=tabulate(influencedNodes(:));
        nodesCount(tab(:,1))=nodesCount(tab(:,1))+tab(:,2);
        for iNode=unique(influencedNodes)'
            nodeConnectElem=sum(QCelemNodeNo==iNode,2);
            cnt=sum(nodeConnectElem);
            nodeConnectElem=find(nodeConnectElem);
            if nodesCount(iNode)>=brokenElemExpandLimit(cnt)
                brokenDOF=[brokenDOF,iNode*2-1,iNode*2];
                if find(QCbrokenLink(nodeConnectElem)==0)
                    forceBroken=true;
                end
                QCbrokenLink(nodeConnectElem)=1;
                
            end
        end
        if forceBroken
            newtoniter=newtoniter-1;
            continue;
        end

        delDOF=unique([delDOF;brokenDOF(:)]);
        brokenDOF=setdiff(delDOF,QCboundaryCfg.dispConstraint);
        saveBrokenValue=disp_init(brokenDOF);
    end


    Kdisp(:,delDOF)=[];
    Kdisp(delDOF,:)=[];
    Rdisp(delDOF)=[];
    %%%%%%%%% Stop crit %%%%%%%%%%
    norm_res=norm(Rdisp);
    if norm_res>10000000
        error('Extreme imbalance');
    end
    if norm_res<computeCfg.critNormRes
        if display
            fprintf('NewtonIter = %2d, Norm of residual = %4.6f, end iteration\n',newtoniter,norm_res);
        end
        break;
    end
    %%%%%%% end stop crit %%%%%%%%
    if ~isempty(brokenElems)
        if strcmp(computeCfg.solveDispMethod,'distributed')
            du=-solve_linear_system(Kdisp,Rdisp,'pcg-dist');
        else
            du=-solve_linear_system(Kdisp,Rdisp,'pcg');
        end
    else
        du=-solve_linear_system(Kdisp,Rdisp,computeCfg.solveDispMethod);
    end
    nanInDu=find(isnan(du));
    infInDu=find(isinf(du));
    if ~isempty(nanInDu) || ~isempty(infInDu)
        disp('(NaN or inf in this step.)');
    end
    du(nanInDu)=0;
    du(infInDu)=0;
    norm_du=norm(du);
    if display && norm_du<1e-10
        fprintf('NewtonIter = %2d, Norm of residual = %4.6f, Norm of du %4.6f, end iteration\n',newtoniter,norm_res,norm_du);
        break;
    end
    % Extremely ill-condition, du might be uncertain for static steps.
    disp_free=full(disp_n);
    disp_free(delDOF)=[];
    disp_free=disp_free+du;
    disp_n=sparse([setdiff(1:QCelemCfg.meshDispDOF,delDOF),QCboundaryCfg.dispConstraint',brokenDOF'],1,...
        [disp_free;saveConsValue;saveBrokenValue],QCelemCfg.meshDispDOF,1);
    %%%%% end Newton-Raphson %%%%%%
    if display
        fprintf('NewtonIter = %2d, Norm of residual = %4.6f, Norm of du %4.6f\n',newtoniter,norm_res,norm_du);
    end
end
if QC_projectCfg.fracture
    fprintf('Max l/Lc %.5f\n',max(lLc));
    QCbrokenLink=tmpBrokenLink;
end
result.disp_net=disp_n;
disp_n1=disp_n;
end
