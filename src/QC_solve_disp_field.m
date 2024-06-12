function [disp_n1,disp_net,result,Eref,canRefine,QCbrokenLink,QCbrokenElem]=QC_solve_disp_field(projectCfg,boundaryCfg,computeCfg,materialCfg,...
    elemCfg,nodeCoord,elemNodeNo,elemNodePos,elemDofLabel,QCelemLength,QCnodeCoord,QCelemNodeNo,QCelemDofLabel,QCelemCfg,QCboundaryCfg,QCbrokenLink,QCbrokenElem,disp_n,N,pNpxi,pNpeta,pNpzeta,M,assem,QCassem,display)
if display
    disp('Solving QC displacement field.');
end
warning('off','MATLAB:singularMatrix');
brokenElemExpandLimit=[1,2,3,4,4,5,6,6];
%%%%%%%%%%%%%%%%%%%%%% Update QC %%%%%%%%%%%%%%%%%%%%%%%%
[netNodeInElem,netElemInElem,elemArea,elemContainNetNode,elemContainNetElem,canRefine]=QC_update_info(elemCfg,QCnodeCoord,QCelemNodeNo,nodeCoord,elemNodeNo,elemNodePos,projectCfg);
QCshapeFuncCoeff=QC_shape_function_coefficient(elemNodePos);
%%%%%%%%%%%%%%%%%%%% End Update QC %%%%%%%%%%%%%%%%%%%%%%
willRefine=zeros(1,elemCfg.nElem);
Eref=[];
lLc=zeros(elemCfg.nElem,1);
tmpBrokenElem=zeros(size(QCbrokenElem));

Voigt=[1,1;2,2;1,2];
meshDispDOF=elemCfg.meshDispDOF;
elemDispDOF=elemCfg.elemDispDOF;
DOF=elemCfg.phyDOF;
nElem=elemCfg.nElem;
nNode=elemCfg.nNode;
nFaceNode=elemCfg.nFaceNode;
nElemNode=elemCfg.nElemNode;
nGaussPoint=elemCfg.nGaussPoint;
W=elemCfg.gaussWeight;
P=computeCfg.P;

kB=materialCfg.kB;
T=materialCfg.T;
b=materialCfg.b;
Lc_L=materialCfg.Lc_L;
Lc=materialCfg.Lc;
% set displacement boundary
[disp_n,saveConsValue]=set_disp_boundary(disp_n,boundaryCfg,computeCfg);

disp_init=disp_n;
newtoniter = 0;
while(newtoniter<=10)

    disp_last=disp_n;
    % Initiation
    %     deformed_nodeCoordRef=nodeCoordRef+disp_n;
    %     elemNodePos=deformed_nodeCoordRef(elemDofLabel);
    newtoniter=newtoniter+1;
    Kdisp_quad=zeros(elemDispDOF^2,nElem);
    Rdisp_quad=zeros(elemDispDOF,nElem);
    disp_net=zeros(length(QCnodeCoord),2);
    elemBrokenLink_quad={};
    %%%%%%%%%%%%% optional %%%%%%%%%%%%%%
    SXX_quad=zeros(nGaussPoint,nElem);
    SYY_quad=zeros(nGaussPoint,nElem);
    SXY_quad=zeros(nGaussPoint,nElem);
    fnet_quad=zeros(length(QCelemNodeNo),1);
    H_quad=zeros(nGaussPoint,nElem);
    Hneg_quad=zeros(nGaussPoint,nElem);
    strainEnergy_quad=zeros(1,nElem);
    %%%%%%%%%%% end optional %%%%%%%%%%%%

    %%%%%%%%%%%%% Interpolate %%%%%%%%%%%%%%
    for iElem=1:nElem
        disp_e=disp_n(elemDofLabel(iElem,:)');
        as=QCshapeFuncCoeff(iElem,1:3);
        bs=QCshapeFuncCoeff(iElem,4:6);
        cs=QCshapeFuncCoeff(iElem,7:9);
        area=elemArea(iElem);
        N1=@(p)(as(1)+bs(1)*p(1)+cs(1)*p(2))/2/area;
        N2=@(p)(as(2)+bs(2)*p(1)+cs(2)*p(2))/2/area;
        N3=@(p)(as(3)+bs(3)*p(1)+cs(3)*p(2))/2/area;

        netNodes=cell2mat(elemContainNetNode(iElem,1));
        for netNodeNo=netNodes'
            %             netNodeNo=netNodes(iNetNode);
            p=QCnodeCoord(netNodeNo,:);
            N1_net=N1(p);    N2_net=N2(p);    N3_net=N3(p);
            Nnet=[N1_net,0,N2_net,0,N3_net,0;0,N1_net,0,N2_net,0,N3_net];
            disp_net(netNodeNo,:)=Nnet*disp_e;
        end
    end
    %%%%%%%%%%% End Interpolate %%%%%%%%%%%%

    %%%%%%%%%%% QC Integrate %%%%%%%%%%%%
    parfor iElem=1:nElem
        if ~(projectCfg.fracture && QCbrokenElem(iElem)==1)
        material=cell2mat(materialCfg.materials(1,1));
        parN=N;
        disp_e=disp_n(elemDofLabel(iElem,:)');

        iGauss=1;
        area=elemArea(iElem);
        N_p=parN(iGauss,:);        
        netElems=cell2mat(elemContainNetElem(iElem,:));
        [J_quad]=elem_jacobian(iElem,elemNodePos,pNpxi,pNpeta,pNpzeta,elemCfg);
        J=reshape(J_quad(:,iGauss),DOF,DOF);
        Gu=compute_Gu(J,iGauss,elemCfg,pNpxi,pNpeta,pNpzeta);
        Farray=Gu*disp_e;
        F=[Farray(1),Farray(3);...
            Farray(4),Farray(2)];
        F=F+eye(DOF,DOF);
        [Bcal,B0,Bu]=compute_B0(F,J,iGauss,elemCfg,pNpxi,pNpeta,pNpzeta);
        detJ=det(J);
        S=zeros(2,2);
        Cvoigt=zeros(3,3);

        elemBrokenLink=zeros(length(netElems),1);
        for iNetElem=1:length(netElems)
            netElemNo=netElems(iNetElem);
            netElemNodeNo=QCelemNodeNo(netElemNo,:);
            netElemDisp=[disp_net(netElemNodeNo(:),:)];
            netElemNodePos=QCnodeCoord(netElemNodeNo,:);
            deformed=netElemNodePos+netElemDisp;
            vec=deformed(1,:)-deformed(2,:);
            origin_vec=netElemNodePos(1,:)-netElemNodePos(2,:);
            l=norm(vec);
            L=QCelemLength(netElemNo);
            Ncurr=vec/norm(vec);
            Norig=origin_vec/norm(origin_vec);
            Lc=L*Lc_L;
            lbd=l/L;
            if newtoniter>1
                lLc(iElem)=max(lLc(iElem),l/Lc);
            end
            if (QCbrokenLink(netNodeNo)==0 && l/Lc < projectCfg.critlbd) || ~projectCfg.fracture
                x=l-L;
%                 x=l;
                if abs(x/Lc-1)<0.01
                    x=0.9*Lc;
                end
                f=kB*T/b*(1/4*(1-x/Lc)^(-2)-1/4+x/Lc);
                k=kB*T/b*(1/2/Lc*(1-x/Lc)^(-3)+1/Lc);
                
                S=S+1/lbd*f*L*kron(Norig',Norig);
                for p=1:3
                    for q=1:3
                        if p>=q
                            Cvoigt(p,q)=Cvoigt(p,q)+(1/lbd^2*k*L^2-1/lbd^3*f*L)*Norig(Voigt(p,1))*Norig(Voigt(p,2))*Norig(Voigt(q,1))*Norig(Voigt(q,2));
                        end
                    end
                end
                
            elseif newtoniter>1
                elemBrokenLink(iNetElem)=1;
            end
        end
        if newtoniter>1  && projectCfg.fracture && length(find(elemBrokenLink))/length(netElems)>projectCfg.critElemFail && ~canRefine(iElem)
            tmpBrokenElem(iElem)=1;
        else
            tmpBrokenElem(iElem)=0;
        end
        elemBrokenLink_quad(iElem)={netElems(find(elemBrokenLink))};

        S=S/area;
        for p=1:3
            for q=1:3
                if p<q
                    Cvoigt(p,q)=Cvoigt(q,p);
                end
            end
        end
        Cvoigt=Cvoigt/area;
        Svoigt=[S(1,1);S(2,2);S(1,2)];

        Kgeo=zeros(DOF*nElemNode,DOF*nElemNode);
        for i=1:nElemNode
            for j=1:nElemNode
                Kij=eye(DOF,DOF);
                Kij=Kij*(W(iGauss)*Bcal(:,i)'*S*Bcal(:,j)*detJ);
                Kgeo((DOF*i-(DOF-1)):(DOF*i),(DOF*j-(DOF-1)):(DOF*j))=Kij;
            end
        end
        Kmat=W(iGauss)*B0'*Cvoigt*B0*detJ;
        R=W(iGauss)*B0'*Svoigt*detJ;

        Kdisp_quad(:,iElem)=Kgeo(:)+Kmat(:);
        Rdisp_quad(:,iElem)=R;
        %%%%%%%%%%%%% optional %%%%%%%%%%%%%%
        %         H_quad(:,iElem)=strainEnergy(1);
        %         Hneg_quad(:,iElem)=strainEnergy(2);
        SXX_quad(:,iElem)=Svoigt(1);
        SYY_quad(:,iElem)=Svoigt(2);
        SXY_quad(:,iElem)=Svoigt(3);
        %         strainEnergy_quad(1,iElem)=W(iGauss)*(s*strainEnergy(1)+strainEnergy(2))*detJ;
        %%%%%%%%%%% end optional %%%%%%%%%%%%
        else
            elemBrokenLink_quad(iElem)={[]};
        end
    end

    %%%%%%%%%%%%%%%%%%%%% Assembling %%%%%%%%%%%%%%%%%%%%%%%%%
    Kdisp=sparse(assem.m8_8i(:),assem.m8_8j(:),Kdisp_quad(:),meshDispDOF,meshDispDOF);
    Rdisp=sparse(assem.m8_1(:),1,Rdisp_quad(:),meshDispDOF,1);
    % Interpolating stress from Gauss points to nodes, the denominator depends
    % on the number of Gauss points for interpolation of each point.
    if nGaussPoint==1
        H_quad=repmat(H_quad,nElemNode,1);
        Hneg_quad=repmat(Hneg_quad,nElemNode,1);
        SXX_quad=repmat(SXX_quad,nElemNode,1);
        SYY_quad=repmat(SYY_quad,nElemNode,1);
        SXY_quad=repmat(SXY_quad,nElemNode,1);
    end

    result.disp_net=disp_net;
    result.SXX=full(sparse(assem.m4_1,1,SXX_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.SYY=full(sparse(assem.m4_1,1,SYY_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.SXY=full(sparse(assem.m4_1,1,SXY_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.fnet=fnet_quad;
    result.strainEnergy=sum(strainEnergy_quad);
    result.reaction=full(Rdisp);
    %%%%%%%%%%%%%%%%% end assembling %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Newton-Raphson %%%%%%%%
    brokenElems=find(QCbrokenElem);
    delDOF=boundaryCfg.dispConstraint;
    brokenDOF=[];
    saveBrokenValue=[];
    
    if projectCfg.fracture && ~isempty(brokenElems)
        forceBroken=false;
        influencedNodes=elemNodeNo(brokenElems,:);
        nodesCount=zeros(elemCfg.nNode,1);
        tab=tabulate(influencedNodes(:));
        nodesCount(tab(:,1))=nodesCount(tab(:,1))+tab(:,2);
        for iNode=unique(influencedNodes)'
            nodeConnectElem=sum(elemNodeNo==iNode,2);
            cnt=sum(nodeConnectElem);
            nodeConnectElem=find(nodeConnectElem);
            if nodesCount(iNode)>=brokenElemExpandLimit(cnt)
                brokenDOF=[brokenDOF,iNode*2-1,iNode*2];
                if find(QCbrokenElem(nodeConnectElem)==0)
                    forceBroken=true;
                end
                QCbrokenElem(nodeConnectElem)=1;
            end
        end
        if forceBroken
            newtoniter=newtoniter-1;
            continue;
        end

        delDOF=unique([delDOF;brokenDOF(:)]);
        brokenDOF=setdiff(delDOF,boundaryCfg.dispConstraint);
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
    disp_n=sparse([setdiff(1:elemCfg.meshDispDOF,delDOF),boundaryCfg.dispConstraint',brokenDOF'],1,...
        [disp_free;saveConsValue;saveBrokenValue],elemCfg.meshDispDOF,1);
    %%%%% end Newton-Raphson %%%%%%
    
    if display
        fprintf('NewtonIter = %2d, Norm of residual = %4.6f, Norm of du %4.6f\n',newtoniter,norm_res,norm_du);
    end
end
%
% parfor iElem=1:QCelemCfg.nElem
%     if ~isempty(find(QCboundaryCfg.dispConstraint==QCelemDofLabel(iElem,1))) || ~isempty(find(QCboundaryCfg.dispConstraint==QCelemDofLabel(iElem,2))) ||...
%        ~isempty(find(QCboundaryCfg.dispConstraint==QCelemDofLabel(iElem,3))) || ~isempty(find(QCboundaryCfg.dispConstraint==QCelemDofLabel(iElem,4)))
%     material=cell2mat(materialCfg.materials(1,1));
% %         Kdisp_e=zeros(elemDispDOF);
%     Rdisp_e=zeros(QCelemCfg.elemDispDOF,1);
%     %%%%%%%%%%%%%%%%%%%%%% QC %%%%%%%%%%%%%%%%%%%%%%%%
%     netElemNodeNo=QCelemNodeNo(iElem,:);
%     netElemDisp=[disp_net(netElemNodeNo(:),:)];
%     netElemNodePos=QCnodeCoord(netElemNodeNo,:);
%     deformed=netElemNodePos+netElemDisp;
%     vec=deformed(1,:)-deformed(2,:);
%     origin_vec=netElemNodePos(1,:)-netElemNodePos(2,:);
%     l=norm(vec);
%     L=QCelemLength(iElem);
%     Ncurr=vec/norm(vec);
%     Norig=origin_vec/norm(origin_vec);
%     lbd=l/L;
%     x=l-L;
%     Lc=L*Lc_L;
%     f=kB*T/b*(1/4*(1-x/Lc)^(-2)-1/4+x/Lc);
%     fnet_quad(iElem,1)=f*b/kB/T;
%     k=kB*T/b*(1/2/Lc*(1-x/Lc)^(-3)+1/Lc);
%
%     S=zeros(2,2);
% %         Cvoigt=zeros(3,3);
%     S=S+1/lbd*f*L*kron(Norig',Norig);
% %         for p=1:3
% %             for q=1:3
% %                 if p>=q
% %                     i=Voigt(p,1);
% %                     j=Voigt(p,2);
% %                     m=Voigt(q,1);
% %                     n=Voigt(q,2);
% %                     Cvoigt(p,q)=Cvoigt(p,q)+(1/lbd^2*k*L^2-1/lbd^3*f*L)*Norig(i)*Norig(j)*Norig(m)*Norig(n);
% %                 end
% %             end
% %         end
% %         for p=1:3
% %             for q=1:3
% %                 if p<q
% %                     Cvoigt(p,q)=Cvoigt(q,p);
% %                 end
% %             end
% %         end
%     Svoigt=[S(1,1);S(2,2);S(1,2)];
%
%     x1=netElemNodePos(1,1); x2=netElemNodePos(2,1);
%     y1=netElemNodePos(1,2); y2=netElemNodePos(2,2);
%
% %         sinT=(y2-y1)/L;
% %         cosT=(x2-x1)/L;
% %         F=lbd*kron(Ncurr',Norig)+0*[sinT^2,-sinT*cosT;-sinT*cosT,cosT^2];
%     F=lbd*kron(Ncurr',Norig);
%     lx=(x2-x1)/2; ly=(y2-y1)/2;
%     pN1px=-lx*2/L^2;
%     pN1py=-ly*2/L^2;
%     pN2px=-pN1px;
%     pN2py=-pN1py;
%     pNpX=[pN1px,pN2px;pN1py,pN2py];
%     Bu=zeros(3,QCelemCfg.elemDispDOF);
%     B0=zeros(3,QCelemCfg.elemDispDOF);
%     for iNode=1:QCelemCfg.nElemNode
%         tmp=[pNpX(1,iNode),0;...
%              0,pNpX(2,iNode);...
%              pNpX(2,iNode),pNpX(1,iNode)];
%         Bu(:,iNode*2-1:iNode*2)=tmp;
%         B0(:,iNode*2-1:iNode*2)=tmp*F';
%     end
%
% %         Bcal=pNpX;
%
%     detJ=1;
%
%     iGauss=1;
% %         Kgeo=zeros(DOF*nElemNode,DOF*nElemNode);
% %         for iC=1:nElemNode
% %             for jC=1:nElemNode
% %                 Kij=eye(DOF,DOF);
% %                 Kij=Kij*(W(iGauss)*Bcal(:,iC)'*S*Bcal(:,jC)*detJ);
% %                 Kgeo((DOF*iC-(DOF-1)):(DOF*iC),(DOF*jC-(DOF-1)):(DOF*jC))=Kij;
% %             end
% %         end
% %         Kdisp_e=Kdisp_e+W(iGauss)*B0'*Cvoigt*B0*detJ+Kgeo;
%     Rdisp_e=Rdisp_e+B0'*Svoigt;
%
% %         Kdisp_quad(:,iElem)=Kdisp_e(:);
%     Rnet_quad(:,iElem)=Rdisp_e;
%     end
% end
% Rnet=sparse(QCassem.m8_1(:),1,Rnet_quad(:),QCelemCfg.meshDispDOF,1);
% result.reaction=full(Rnet);
if  projectCfg.fracture
fprintf('Max l/Lc %.5f\n',max(lLc));

brokenLink=[];
for iElem=1:nElem
    tmp=cell2mat(elemBrokenLink_quad(iElem));
    brokenLink=[brokenLink;tmp];
    if ~isempty(tmp) && canRefine(iElem)
        willRefine(iElem)=1;
    end
end
Eref=find(willRefine);
QCbrokenLink(unique(brokenLink))=1;
QCbrokenElem=max(QCbrokenElem,tmpBrokenElem);
tmp=elemContainNetElem(find(QCbrokenElem),1);
netElemInBrokenElem=[];
for itmp=1:length(tmp)
    netElemInBrokenElem=[netElemInBrokenElem;cell2mat(tmp(itmp))];
end
QCbrokenLink(netElemInBrokenElem(:))=1;
end
disp_n1=disp_n;
end

