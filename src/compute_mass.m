function [M]=compute_mass(elemCfg,elemNodePos,N,pNpxi,pNpeta,pNpzeta,materialCfg,assem)
% Compute mass matrix which is diagnose.

% The mass matrix is diag. For simplification, it is treated as a column.
% Me=int(rho*Ndisp^T)
% In the loop, M is only (nNode x 1) for nodes, but the masses of the node
% for both directions are the same, so there is function kron for duplicate
% the mass for each node.
disp('Compute mass matrix')
%%%%%%%%%%% Compare parfor with unpolished code %%%%%%%%%%%%
% M=sparse(elemCfg.nNode,1);
% for iElem=1:elemCfg.nElem
%     Me=zeros(4,1);
%     for iGauss=1:elemCfg.nGaussPoint
%         J=Jacobian(iElem,iGauss);
%         Me=Me+rho*Ndisp(iGauss,:)'*det(J);
%     end
%     M(elemNodeNo(iElem,:))= M(elemNodeNo(iElem,:))+Me(:);
% end
% M=kron(M,ones(elemCfg.phyDOF,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mquad=zeros(elemCfg.nElemNode,elemCfg.nElem);
nElem=elemCfg.nElem;
nGaussPoint=elemCfg.nGaussPoint;
nElemNode=elemCfg.nElemNode;
DOF=elemCfg.phyDOF;
weight=elemCfg.gaussWeight;
parfor iElem=1:nElem
    material=cell2mat(materialCfg.materials(materialCfg.elemMatType(iElem),1));
    rho=material.rho;

    parN=N; %reduce communication cost
    Me=zeros(nElemNode,1);
    [J_quad]=elem_jacobian(iElem,elemNodePos,pNpxi,pNpeta,pNpzeta,elemCfg);
    for iGauss=1:nGaussPoint
        J=reshape(J_quad(:,iGauss),DOF,DOF);
        Me=Me+weight(iGauss)*rho*parN(iGauss,:)'*det(J);
    end
    Mquad(:,iElem)=Me;
end
%%%%%%%%%%%%% unpolished assembling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for iElem=1:nElem
%     M(elemNodeNo(iElem,:))= M(elemNodeNo(iElem,:))+Mquad(:,iElem);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=sparse(assem.m4_1,1,Mquad,elemCfg.nNode,1);
M=kron(M,ones(elemCfg.phyDOF,1));