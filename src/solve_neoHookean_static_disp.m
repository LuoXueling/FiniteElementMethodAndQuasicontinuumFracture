function [disp_n1,velo_n1,acce_n1,result]=solve_neoHookean_static_disp(projectCfg,boundaryCfg,computeCfg,materialCfg,...
elemCfg,elemNodeNo,elemNodePos,elemDofLabel,disp_n,velo_n,acce_n,N,pNpxi,pNpeta,pNpzeta,M,assem,display)

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

[f_face_quad,result]=integrate_face_force(projectCfg,boundaryCfg,computeCfg,elemCfg,elemNodePos);

if display
disp(cell2mat(strcat('Solving neo-Hookean static-disp with method:',{32},computeCfg.solveDispMethod)));
end
%     du=sparse(elemCfg.meshDispDOF-length(boundaryCfg.dispConstraint),1);
% set displacement boundary
[disp_n,saveConsValue]=set_disp_boundary(disp_n,boundaryCfg,computeCfg);
disp_init=disp_n;
% Save constrained displacement to recover after each iteration
iter=0;
while(iter<=50)
    % Initiation
    iter=iter+1;
    Kdisp_quad=zeros(elemDispDOF^2,nElem);
    Rdisp_quad=zeros(elemDispDOF,nElem);
    %%%%%%%%%%%%% optional %%%%%%%%%%%%%%
    sigmaXX_quad=zeros(nGaussPoint,nElem);
    sigmaYY_quad=zeros(nGaussPoint,nElem);
    sigmaXY_quad=zeros(nGaussPoint,nElem);
    epsilonXX_quad=zeros(nGaussPoint,nElem);
    epsilonYY_quad=zeros(nGaussPoint,nElem);
    epsilonXY_quad=zeros(nGaussPoint,nElem);
    SXX_quad=zeros(nGaussPoint,nElem);
    SYY_quad=zeros(nGaussPoint,nElem);
    SXY_quad=zeros(nGaussPoint,nElem);
    PXX_quad=zeros(nGaussPoint,nElem);
    PYY_quad=zeros(nGaussPoint,nElem);
    PXY_quad=zeros(nGaussPoint,nElem);
    if DOF==3
        sigmaZZ_quad=zeros(nGaussPoint,nElem);
        sigmaXZ_quad=zeros(nGaussPoint,nElem);
        sigmaYZ_quad=zeros(nGaussPoint,nElem);
        epsilonZZ_quad=zeros(nGaussPoint,nElem);
        epsilonXZ_quad=zeros(nGaussPoint,nElem);
        epsilonYZ_quad=zeros(nGaussPoint,nElem);
        SZZ_quad=zeros(nGaussPoint,nElem);
        SXZ_quad=zeros(nGaussPoint,nElem);
        SYZ_quad=zeros(nGaussPoint,nElem);
        PZZ_quad=zeros(nGaussPoint,nElem);
        PXZ_quad=zeros(nGaussPoint,nElem);
        PYZ_quad=zeros(nGaussPoint,nElem);
    end
    strainEnergy_quad=zeros(1,nElem);
    %%%%%%%%%%% end optional %%%%%%%%%%%%
    parfor iElem=1:nElem
        material=cell2mat(materialCfg.materials(materialCfg.elemMatType(iElem),1));
        
        parN=N;
        disp_e=disp_n(elemDofLabel(iElem,:)');
        [J_quad]=elem_jacobian(iElem,elemNodePos,pNpxi,pNpeta,pNpzeta,elemCfg);
        Kdisp_e=zeros(elemDispDOF);
        Rdisp_e=zeros(elemDispDOF,1);
        %%%%%%%%%%%%% optional %%%%%%%%%%%%%%
        sigmaXX_e=zeros(nGaussPoint,1);
        sigmaYY_e=zeros(nGaussPoint,1);
        sigmaXY_e=zeros(nGaussPoint,1);
        sigmaZZ_e=zeros(nGaussPoint,1);
        sigmaXZ_e=zeros(nGaussPoint,1);
        sigmaYZ_e=zeros(nGaussPoint,1);
        SXX_e=zeros(nGaussPoint,1);
        SYY_e=zeros(nGaussPoint,1);
        SXY_e=zeros(nGaussPoint,1);
        SZZ_e=zeros(nGaussPoint,1);
        SXZ_e=zeros(nGaussPoint,1);
        SYZ_e=zeros(nGaussPoint,1);
        PXX_e=zeros(nGaussPoint,1);
        PYY_e=zeros(nGaussPoint,1);
        PXY_e=zeros(nGaussPoint,1);
        PZZ_e=zeros(nGaussPoint,1);
        PXZ_e=zeros(nGaussPoint,1);
        PYZ_e=zeros(nGaussPoint,1);
        epsilonXX_e=zeros(nGaussPoint,1);
        epsilonYY_e=zeros(nGaussPoint,1);
        epsilonXY_e=zeros(nGaussPoint,1);
        epsilonZZ_e=zeros(nGaussPoint,1);
        epsilonXZ_e=zeros(nGaussPoint,1);
        epsilonYZ_e=zeros(nGaussPoint,1);
        strainEnergy_e=0;
        %%%%%%%%%%% end optional %%%%%%%%%%%%
        for iGauss=1:nGaussPoint
            J=reshape(J_quad(:,iGauss),DOF,DOF);
            
            Gu=compute_Gu(J,iGauss,elemCfg,pNpxi,pNpeta,pNpzeta);
            Farray=Gu*disp_e;
            F=zeros(DOF,DOF);
            if DOF==2
                F=[Farray(1),Farray(3);...
                   Farray(4),Farray(2)];
            elseif DOF==3
                F=[Farray(1),Farray(4),Farray(7);...
                   Farray(8),Farray(2),Farray(5);...
                   Farray(6),Farray(9),Farray(3)];
            end
            F=F+eye(DOF,DOF);
            
            [Bcal,B0,Bu]=compute_B0(F,J,iGauss,elemCfg,pNpxi,pNpeta,pNpzeta);
            
            detJ=det(J);

            s=1;
            [S,Svoigt,Pvoigt,sigma,epsilon,strainEnergy,matStiffDensity] =...
                elastic_neoHookean(F,s,material);
            %%%%%%%%%%%%% Geo matrix %%%%%%%%%%%%%%
            Kgeo=zeros(DOF*nElemNode,DOF*nElemNode);
            for i=1:nElemNode
                for j=1:nElemNode
                    Kij=eye(DOF,DOF);
                    Kij=Kij*(W(iGauss)*Bcal(:,i)'*S*Bcal(:,j)*detJ);
                    Kgeo((DOF*i-(DOF-1)):(DOF*i),(DOF*j-(DOF-1)):(DOF*j))=Kij;
                end
            end
            %%%%%%%%%%% End Geo matrix %%%%%%%%%%%%
            %%%%%%%%%%%%% Mat matrix %%%%%%%%%%%%%%
            Kmat=W(iGauss)*B0'*matStiffDensity*B0*detJ;
            %%%%%%%%%%% End Mat matrix %%%%%%%%%%%%
            % Stiffness matrix and residual vector
            Kdisp_e=Kdisp_e+Kgeo+Kmat;
            Rdisp_e=Rdisp_e+W(iGauss)*B0'*Svoigt*detJ;
            %%%%%%%%%%%%% optional %%%%%%%%%%%%%%
            strainEnergy_e=strainEnergy_e+W(iGauss)*(s*strainEnergy(1)+strainEnergy(2))*detJ;
            if DOF==2
                sigmaXX_e(iGauss)=sigma(1);
                sigmaYY_e(iGauss)=sigma(2);
                sigmaXY_e(iGauss)=sigma(3);
                SXX_e(iGauss)=Svoigt(1);
                SYY_e(iGauss)=Svoigt(2);
                SXY_e(iGauss)=Svoigt(3);
                PXX_e(iGauss)=Pvoigt(1);
                PYY_e(iGauss)=Pvoigt(2);
                PXY_e(iGauss)=Pvoigt(3);
                epsilonXX_e(iGauss)=epsilon(1);
                epsilonYY_e(iGauss)=epsilon(2);
                epsilonXY_e(iGauss)=epsilon(3);
            elseif DOF==3
                sigmaXX_e(iGauss)=sigma(1);
                sigmaYY_e(iGauss)=sigma(2);
                sigmaZZ_e(iGauss)=sigma(3);
                sigmaXY_e(iGauss)=sigma(4);
                sigmaYZ_e(iGauss)=sigma(5);
                sigmaXZ_e(iGauss)=sigma(6);
                SXX_e(iGauss)=Svoigt(1);
                SYY_e(iGauss)=Svoigt(2);
                SZZ_e(iGauss)=Svoigt(3);
                SXY_e(iGauss)=Svoigt(4);
                SYZ_e(iGauss)=Svoigt(5);
                SXZ_e(iGauss)=Svoigt(6);
                PXX_e(iGauss)=Pvoigt(1);
                PYY_e(iGauss)=Pvoigt(2);
                PZZ_e(iGauss)=Pvoigt(3);
                PXY_e(iGauss)=Pvoigt(4);
                PYZ_e(iGauss)=Pvoigt(5);
                PXZ_e(iGauss)=Pvoigt(6);
                epsilonXX_e(iGauss)=epsilon(1);
                epsilonYY_e(iGauss)=epsilon(2);
                epsilonZZ_e(iGauss)=epsilon(3);
                epsilonXY_e(iGauss)=epsilon(4);
                epsilonYZ_e(iGauss)=epsilon(5);
                epsilonXZ_e(iGauss)=epsilon(6);
            end
            %%%%%%%%%%% end optional %%%%%%%%%%%%
        end
        Kdisp_quad(:,iElem)=Kdisp_e(:);
        Rdisp_quad(:,iElem)=Rdisp_e;
        %%%%%%%%%%%%% optional %%%%%%%%%%%%%%
        sigmaXX_quad(:,iElem)=sigmaXX_e;
        sigmaYY_quad(:,iElem)=sigmaYY_e;
        sigmaXY_quad(:,iElem)=sigmaXY_e;
        SXX_quad(:,iElem)=SXX_e;
        SYY_quad(:,iElem)=SYY_e;
        SXY_quad(:,iElem)=SXY_e;
        PXX_quad(:,iElem)=PXX_e;
        PYY_quad(:,iElem)=PYY_e;
        PXY_quad(:,iElem)=PXY_e;
        epsilonXX_quad(:,iElem)=epsilonXX_e;
        epsilonYY_quad(:,iElem)=epsilonYY_e;
        epsilonXY_quad(:,iElem)=epsilonXY_e;
        if DOF==3
            sigmaZZ_quad(:,iElem)=sigmaZZ_e;
            sigmaYZ_quad(:,iElem)=sigmaYZ_e;
            sigmaXZ_quad(:,iElem)=sigmaXZ_e;
            SZZ_quad(:,iElem)=SZZ_e;
            SYZ_quad(:,iElem)=SYZ_e;
            SXZ_quad(:,iElem)=SXZ_e;
            PZZ_quad(:,iElem)=PZZ_e;
            PYZ_quad(:,iElem)=PYZ_e;
            PXZ_quad(:,iElem)=PXZ_e;
            epsilonZZ_quad(:,iElem)=epsilonZZ_e;
            epsilonYZ_quad(:,iElem)=epsilonYZ_e;
            epsilonXZ_quad(:,iElem)=epsilonXZ_e;
        end
        strainEnergy_quad(1,iElem)=strainEnergy_e;
        %%%%%%%%%%% end optional %%%%%%%%%%%%
    end
    if projectCfg.surfaceTraction
        Rdisp_quad=Rdisp_quad-f_face_quad;
    end
    %%%%%%%%%%%%%%%%%%%%% Assembling %%%%%%%%%%%%%%%%%%%%%%%%%
    Kdisp=sparse(assem.m8_8i(:),assem.m8_8j(:),Kdisp_quad(:),meshDispDOF,meshDispDOF);
    Rdisp=sparse(assem.m8_1(:),1,Rdisp_quad(:),meshDispDOF,1);
    % Interpolating stress from Gauss points to nodes, the denominator depends
    % on the number of Gauss points for interpolation of each point.
    if nGaussPoint==1
        sigmaXX_quad=repmat(sigmaXX_quad,nElemNode,1);
        sigmaYY_quad=repmat(sigmaYY_quad,nElemNode,1);
        sigmaXY_quad=repmat(sigmaXY_quad,nElemNode,1);
        SXX_quad=repmat(SXX_quad,nElemNode,1);
        SYY_quad=repmat(SYY_quad,nElemNode,1);
        SXY_quad=repmat(SXY_quad,nElemNode,1);
        PXX_quad=repmat(PXX_quad,nElemNode,1);
        PYY_quad=repmat(PYY_quad,nElemNode,1);
        PXY_quad=repmat(PXY_quad,nElemNode,1);
        epsilonXX_quad=repmat(epsilonXX_quad,nElemNode,1);
        epsilonYY_quad=repmat(epsilonYY_quad,nElemNode,1);
        epsilonXY_quad=repmat(epsilonXY_quad,nElemNode,1);
        if DOF==3
            sigmaZZ_quad=repmat(sigmaZZ_quad,nElemNode,1);
            sigmaXZ_quad=repmat(sigmaXZ_quad,nElemNode,1);
            sigmaYZ_quad=repmat(sigmaYZ_quad,nElemNode,1);
            SZZ_quad=repmat(SZZ_quad,nElemNode,1);
            SXZ_quad=repmat(SXZ_quad,nElemNode,1);
            SYZ_quad=repmat(SYZ_quad,nElemNode,1);
            PZZ_quad=repmat(PZZ_quad,nElemNode,1);
            PXZ_quad=repmat(PXZ_quad,nElemNode,1);
            PYZ_quad=repmat(PYZ_quad,nElemNode,1);
            epsilonZZ_quad=repmat(epsilonZZ_quad,nElemNode,1);
            epsilonXZ_quad=repmat(epsilonXZ_quad,nElemNode,1);
            epsilonYZ_quad=repmat(epsilonYZ_quad,nElemNode,1);
        end
    end
    result.sigmaXX=full(sparse(assem.m4_1,1,sigmaXX_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.sigmaYY=full(sparse(assem.m4_1,1,sigmaYY_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.sigmaXY=full(sparse(assem.m4_1,1,sigmaXY_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.sigmaXX_quad=sigmaXX_quad;
    result.sigmaYY_quad=sigmaYY_quad;
    result.sigmaXY_quad=sigmaXY_quad;
    result.SXX=full(sparse(assem.m4_1,1,SXX_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.SYY=full(sparse(assem.m4_1,1,SYY_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.SXY=full(sparse(assem.m4_1,1,SXY_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.PXX=full(sparse(assem.m4_1,1,PXX_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.PYY=full(sparse(assem.m4_1,1,PYY_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.PXY=full(sparse(assem.m4_1,1,PXY_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.epsilonXX=full(sparse(assem.m4_1,1,epsilonXX_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.epsilonYY=full(sparse(assem.m4_1,1,epsilonYY_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.epsilonXY=full(sparse(assem.m4_1,1,epsilonXY_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    if DOF==3
        result.sigmaZZ=full(sparse(assem.m4_1,1,sigmaZZ_quad(:),nNode,1))...
            ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
        result.sigmaYZ=full(sparse(assem.m4_1,1,sigmaYZ_quad(:),nNode,1))...
            ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
        result.sigmaXZ=full(sparse(assem.m4_1,1,sigmaXZ_quad(:),nNode,1))...
            ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
        result.sigmaZZ_quad=sigmaZZ_quad;
        result.sigmaYZ_quad=sigmaYZ_quad;
        result.sigmaXZ_quad=sigmaXZ_quad;
        result.SZZ=full(sparse(assem.m4_1,1,SZZ_quad(:),nNode,1))...
            ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
        result.SYZ=full(sparse(assem.m4_1,1,SYZ_quad(:),nNode,1))...
            ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
        result.SXZ=full(sparse(assem.m4_1,1,SXZ_quad(:),nNode,1))...
            ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
        result.PZZ=full(sparse(assem.m4_1,1,PZZ_quad(:),nNode,1))...
            ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
        result.PYZ=full(sparse(assem.m4_1,1,PYZ_quad(:),nNode,1))...
            ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
        result.PXZ=full(sparse(assem.m4_1,1,PXZ_quad(:),nNode,1))...
            ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
        result.epsilonZZ=full(sparse(assem.m4_1,1,epsilonZZ_quad(:),nNode,1))...
            ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
        result.epsilonYZ=full(sparse(assem.m4_1,1,epsilonYZ_quad(:),nNode,1))...
            ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
        result.epsilonXZ=full(sparse(assem.m4_1,1,epsilonXZ_quad(:),nNode,1))...
            ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    end

    result.strainEnergy=sum(strainEnergy_quad);
    result.reaction=full(Rdisp);
    %%%%%%%%%%%%%%%%% end assembling %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Newton-Raphson %%%%%%%%
    Kdisp(:,boundaryCfg.dispConstraint)=[];
    Kdisp(boundaryCfg.dispConstraint,:)=[];
    Rdisp(boundaryCfg.dispConstraint(:))=[];
    %%%%%%%%% Stop crit %%%%%%%%%%
    norm_res=norm(Rdisp);
    if norm_res<computeCfg.critNormRes
        if display
        fprintf('Step = %2d, Norm of residual = %4.6f, end iteration\n',iter,norm_res);
        end
        break;
    end
    if norm_res>10000000 && iter>20
        error('Extreme imbalance.');
    end
    %%%%%%% end stop crit %%%%%%%%
    du=-solve_linear_system(Kdisp,Rdisp,computeCfg.solveDispMethod);
    % Extremely ill-condition, du might be uncertain for static steps.
    disp_free=full(disp_n);
    disp_free(boundaryCfg.dispConstraint)=[];
    disp_free=disp_free+1*du;
    disp_n=sparse([boundaryCfg.dispFree,boundaryCfg.dispConstraint'],1,...
        [disp_free;saveConsValue],elemCfg.meshDispDOF,1);
    %%%%% end Newton-Raphson %%%%%%
    norm_du=norm(du);
    if display
    fprintf('Step = %2d, Norm of residual = %4.6f, Norm of du %4.6f\n',iter,norm_res,norm_du);
    end
end
disp_n1=disp_n;
velo_n1=velo_n;
acce_n1=acce_n;
result.dDisp=disp_n-disp_init;
end
