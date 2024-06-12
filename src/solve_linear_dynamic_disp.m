function [disp_n1,velo_n1,acce_n1,result]=solve_linear_dynamic_disp(projectCfg,boundaryCfg,computeCfg,materialCfg,...
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
disp('Solving linear elastic dynamic-disp');
end
% velocity change for half time step 
velo_half=velo_n+0.5*acce_n*computeCfg.dynamicDt;
% enforce velocity boundary
velo_half=set_velo_boundary(velo_half,boundaryCfg,computeCfg.dynamicVelo);
% displacement change for full time step
disp_n1=disp_n+computeCfg.dynamicDt*velo_half;
% force -> Accelaration
% The following integral is almost the same as the static part,
% no degradation and no stiffness matrix
f_quad=zeros(elemDispDOF,nElem);
%%%%%%%%%%%%% optional %%%%%%%%%%%%%%
sigmaXX_quad=zeros(nGaussPoint,nElem);
sigmaYY_quad=zeros(nGaussPoint,nElem);
sigmaXY_quad=zeros(nGaussPoint,nElem);
epsilonXX_quad=zeros(nGaussPoint,nElem);
epsilonYY_quad=zeros(nGaussPoint,nElem);
epsilonXY_quad=zeros(nGaussPoint,nElem);
if DOF==3
    sigmaZZ_quad=zeros(nGaussPoint,nElem);
    sigmaXZ_quad=zeros(nGaussPoint,nElem);
    sigmaYZ_quad=zeros(nGaussPoint,nElem);
    epsilonZZ_quad=zeros(nGaussPoint,nElem);
    epsilonXZ_quad=zeros(nGaussPoint,nElem);
    epsilonYZ_quad=zeros(nGaussPoint,nElem);
end
strainEnergy_quad=zeros(1,nElem);
%%%%%%%%%%% end optional %%%%%%%%%%%%
parfor iElem=1:nElem
    material=cell2mat(materialCfg.materials(materialCfg.elemMatType(iElem),1));
    rho=material.rho;
    parN=N;
    acce_e=acce_n(elemDofLabel(iElem,:)');
    disp_e=disp_n1(elemDofLabel(iElem,:)'); %%%%%%disp_n1 here
    [J_quad]=elem_jacobian(iElem,elemNodePos,pNpxi,pNpeta,pNpzeta,elemCfg);
    f_e=zeros(elemDispDOF,1);
    %%%%%%%%%%%%% optional %%%%%%%%%%%%%%
    sigmaXX_e=zeros(nGaussPoint,1);
    sigmaYY_e=zeros(nGaussPoint,1);
    sigmaXY_e=zeros(nGaussPoint,1);
    sigmaZZ_e=zeros(nGaussPoint,1);
    sigmaXZ_e=zeros(nGaussPoint,1);
    sigmaYZ_e=zeros(nGaussPoint,1);
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
        Bu=compute_Bu(J,iGauss,elemCfg,pNpxi,pNpeta,pNpzeta);
        Gu=compute_Gu(J,iGauss,elemCfg,pNpxi,pNpeta,pNpzeta);
        detJ=det(J);
        epsilon=Bu*disp_e;
        Farray=Gu*disp_e;
        s=1;
        [sigma,strainEnergy,~] =...
                elastic_linear(epsilon,s,material);
        %%%%%%%%%%%%% optional %%%%%%%%%%%%%%
        strainEnergy_e=strainEnergy_e+W(iGauss)*(s*strainEnergy(1)+strainEnergy(2))*detJ;
        if DOF==2
            sigmaXX_e(iGauss)=sigma(1);
            sigmaYY_e(iGauss)=sigma(2);
            sigmaXY_e(iGauss)=sigma(3);
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
            epsilonXX_e(iGauss)=epsilon(1);
            epsilonYY_e(iGauss)=epsilon(2);
            epsilonZZ_e(iGauss)=epsilon(3);
            epsilonXY_e(iGauss)=epsilon(4);
            epsilonYZ_e(iGauss)=epsilon(5);
            epsilonXZ_e(iGauss)=epsilon(6);
        end
        %%%%%%%%%%% end optional %%%%%%%%%%%%
        % Internal vector
        f_e=f_e+W(iGauss)*Bu'*sigma*detJ;
    end
    f_quad(:,iElem)=f_e;
    %%%%%%%%%%%%% optional %%%%%%%%%%%%%%
    sigmaXX_quad(:,iElem)=sigmaXX_e;
    sigmaYY_quad(:,iElem)=sigmaYY_e;
    sigmaXY_quad(:,iElem)=sigmaXY_e;
    epsilonXX_quad(:,iElem)=epsilonXX_e;
    epsilonYY_quad(:,iElem)=epsilonYY_e;
    epsilonXY_quad(:,iElem)=epsilonXY_e;
    if DOF==3
        sigmaZZ_quad(:,iElem)=sigmaZZ_e;
        sigmaYZ_quad(:,iElem)=sigmaYZ_e;
        sigmaXZ_quad(:,iElem)=sigmaXZ_e;
        epsilonZZ_quad(:,iElem)=epsilonZZ_e;
        epsilonYZ_quad(:,iElem)=epsilonYZ_e;
        epsilonXZ_quad(:,iElem)=epsilonXZ_e;
    end
    strainEnergy_quad(1,iElem)=strainEnergy_e;
    %%%%%%%%%%% end optional %%%%%%%%%%%%
end
if projectCfg.surfaceTraction
    f_quad=f_quad-f_face_quad;
end
%%%%%%%%%%%%%%%%%%%%% Assembling %%%%%%%%%%%%%%%%%%%%%%%%%
f=real(sparse(assem.m8_1,1,f_quad(:),meshDispDOF,1));
result.reaction=full(f);
if nGaussPoint==1
    sigmaXX_quad=repmat(sigmaXX_quad,nElemNode,1);
    sigmaYY_quad=repmat(sigmaYY_quad,nElemNode,1);
    sigmaXY_quad=repmat(sigmaXY_quad,nElemNode,1);
    epsilonXX_quad=repmat(epsilonXX_quad,nElemNode,1);
    epsilonYY_quad=repmat(epsilonYY_quad,nElemNode,1);
    epsilonXY_quad=repmat(epsilonXY_quad,nElemNode,1);
    if DOF==3
        sigmaZZ_quad=repmat(sigmaZZ_quad,nElemNode,1);
        sigmaXZ_quad=repmat(sigmaXZ_quad,nElemNode,1);
        sigmaYZ_quad=repmat(sigmaYZ_quad,nElemNode,1);
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
    result.sigmaXX_quad=sigmaXX_quad;
    result.sigmaYY_quad=sigmaYY_quad;
    result.sigmaXY_quad=sigmaXY_quad;
    result.epsilonZZ=full(sparse(assem.m4_1,1,epsilonZZ_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.epsilonYZ=full(sparse(assem.m4_1,1,epsilonYZ_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.epsilonXZ=full(sparse(assem.m4_1,1,epsilonXZ_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
end
result.strainEnergy=sum(strainEnergy_quad);
%%%%%%%%%%%%%%%%% end assembling %%%%%%%%%%%%%%%%%%%%%%%%
% Cause the Mass matrix is diagnose, simply divide ./ 
acce_n1=-f./M;
acce_n1=set_acce_boundary(acce_n1,boundaryCfg);
% velocity change for the other half time step 
velo_n1=velo_half+0.5*acce_n1*computeCfg.dynamicDt;
% enforce velocity and acceleration boundary
velo_n1=set_velo_boundary(velo_n1,boundaryCfg,computeCfg.dynamicVelo);
end

