function [disp_n1,velo_n1,acce_n1,result]=solve_linear_static_disp(projectCfg,boundaryCfg,computeCfg,materialCfg,...
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
disp(cell2mat(strcat('Solving linear elastic static-disp with method:',{32},computeCfg.solveDispMethod)));
end
%     du=sparse(elemCfg.meshDispDOF-length(boundaryCfg.dispConstraint),1);
% set displacement boundary
[disp_n,saveConsValue]=set_disp_boundary(disp_n,boundaryCfg,computeCfg);
disp_init=disp_n;
% Save constrained displacement to recover after each iteration
iter=0;
while(iter<=20)
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
    FXX_quad=zeros(nGaussPoint,nElem);
    FYX_quad=zeros(nGaussPoint,nElem);
    if DOF==3
        sigmaZZ_quad=zeros(nGaussPoint,nElem);
        sigmaXZ_quad=zeros(nGaussPoint,nElem);
        sigmaYZ_quad=zeros(nGaussPoint,nElem);
        epsilonZZ_quad=zeros(nGaussPoint,nElem);
        epsilonXZ_quad=zeros(nGaussPoint,nElem);
        epsilonYZ_quad=zeros(nGaussPoint,nElem);
    end

    strainEnergyDensity_quad=zeros(nGaussPoint,nElem);
    strainEnergy_quad=zeros(1,nElem);
    stress_quad=zeros(nGaussPoint,nElem);
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
        epsilonXX_e=zeros(nGaussPoint,1);
        epsilonYY_e=zeros(nGaussPoint,1);
        epsilonXY_e=zeros(nGaussPoint,1);
        epsilonZZ_e=zeros(nGaussPoint,1);
        epsilonXZ_e=zeros(nGaussPoint,1);
        epsilonYZ_e=zeros(nGaussPoint,1);
        FXX_e=zeros(nGaussPoint,1);
        FYX_e=zeros(nGaussPoint,1);
        strainEnergyDensity_e=zeros(nGaussPoint,1);
        strainEnergy_e=0;
        stress_e=zeros(nGaussPoint,1);
        %%%%%%%%%%% end optional %%%%%%%%%%%%
        for iGauss=1:nGaussPoint
            J=reshape(J_quad(:,iGauss),DOF,DOF);
            Bu=compute_Bu(J,iGauss,elemCfg,pNpxi,pNpeta,pNpzeta);
            Gu=compute_Gu(J,iGauss,elemCfg,pNpxi,pNpeta,pNpzeta);
            detJ=det(J);
            epsilon=Bu*disp_e;
            Farray=Gu*disp_e;
            s=1;
            [sigma,strainEnergy,stiffnessDensity] =...
                elastic_linear(epsilon,s,material);
            %%%%%%%%%%%%% optional %%%%%%%%%%%%%%
            strainEnergyDensity_e(iGauss)=s*strainEnergy(1)+strainEnergy(2);
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
            FXX_e(iGauss)=Farray(1);
            FYX_e(iGauss)=Farray(4);
            stress_e(iGauss)=(sigma(1)+sigma(2))/2;
            %%%%%%%%%%% end optional %%%%%%%%%%%%
            % Stiffness matrix and residual vector
            Kdisp_e=Kdisp_e+W(iGauss)*Bu'*stiffnessDensity*Bu*detJ;
            Rdisp_e=Rdisp_e+W(iGauss)*Bu'*sigma*detJ;
        end
        Kdisp_quad(:,iElem)=Kdisp_e(:);
        Rdisp_quad(:,iElem)=Rdisp_e;
        %%%%%%%%%%%%% optional %%%%%%%%%%%%%%
        sigmaXX_quad(:,iElem)=sigmaXX_e;
        sigmaYY_quad(:,iElem)=sigmaYY_e;
        sigmaXY_quad(:,iElem)=sigmaXY_e;
        epsilonXX_quad(:,iElem)=epsilonXX_e;
        epsilonYY_quad(:,iElem)=epsilonYY_e;
        epsilonXY_quad(:,iElem)=epsilonXY_e;
        FXX_quad(:,iElem)=FXX_e;
        FYX_quad(:,iElem)=FYX_e;
        if DOF==3
            sigmaZZ_quad(:,iElem)=sigmaZZ_e;
            sigmaYZ_quad(:,iElem)=sigmaYZ_e;
            sigmaXZ_quad(:,iElem)=sigmaXZ_e;
            epsilonZZ_quad(:,iElem)=epsilonZZ_e;
            epsilonYZ_quad(:,iElem)=epsilonYZ_e;
            epsilonXZ_quad(:,iElem)=epsilonXZ_e;
        end
        strainEnergyDensity_quad(:,iElem)=strainEnergyDensity_e;
        strainEnergy_quad(iElem)=strainEnergy_e;
        stress_quad(:,iElem)=stress_e;
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
        epsilonXX_quad=repmat(epsilonXX_quad,nElemNode,1);
        epsilonYY_quad=repmat(epsilonYY_quad,nElemNode,1);
        epsilonXY_quad=repmat(epsilonXY_quad,nElemNode,1);
        FXX_quad=repmat(FXX_quad,nElemNode,1);
        FYX_quad=repmat(FYX_quad,nElemNode,1);
        stress_quad=repmat(stress_quad,nElemNode,1);
        strainEnergyDensity_quad=repmat(strainEnergyDensity_quad,nElemNode,1);
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
    result.FXX=full(sparse(assem.m4_1,1,FXX_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.FYX=full(sparse(assem.m4_1,1,FYX_quad(:),nNode,1))...
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
        result.epsilonZZ=full(sparse(assem.m4_1,1,epsilonZZ_quad(:),nNode,1))...
            ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
        result.epsilonYZ=full(sparse(assem.m4_1,1,epsilonYZ_quad(:),nNode,1))...
            ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
        result.epsilonXZ=full(sparse(assem.m4_1,1,epsilonXZ_quad(:),nNode,1))...
            ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    end
    if DOF==2
        result.hydrostaticStress=(result.sigmaXX+result.sigmaYY)/2;
    else
        result.hydrostaticStress=(result.sigmaXX+result.sigmaYY+result.sigmaZZ)/3;
    end
    result.stress=full(sparse(assem.m4_1,1,stress_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
    result.strainEnergyDensity=full(sparse(assem.m4_1,1,strainEnergyDensity_quad(:),nNode,1))...
        ./sparse(elemNodeNo',1,ones(nElemNode*nElem,1),nNode,1);
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
    if norm_res>1000000 && iter>5
        error('Divergence');
    end
    %%%%%%% end stop crit %%%%%%%%
    du=-solve_linear_system(Kdisp,Rdisp,computeCfg.solveDispMethod);
    % Extremely ill-condition, du might be uncertain for static steps.
    disp_free=full(disp_n);
    disp_free(boundaryCfg.dispConstraint)=[];
    disp_free=disp_free+du;
    disp_n=sparse([boundaryCfg.dispFree,boundaryCfg.dispConstraint'],1,...
        [disp_free;saveConsValue],elemCfg.meshDispDOF,1);
    %%%%% end Newton-Raphson %%%%%%
    if display
    fprintf('Step = %2d, Norm of residual = %4.6f, Norm of du = %4.6f.\n',iter,norm_res,norm(du));
    end
end
disp_n1=disp_n;
velo_n1=velo_n;
acce_n1=acce_n;
result.dDisp=disp_n-disp_init;
end

