function [J_quad]=elem_jacobian(iElem,elemNodePos,pNpxi,pNpeta,pNpzeta,elemCfg)
% Return Jacobians for all Gauss points of the current element. 
% Jacobian will be computed on every gauss point, so using Jacobian in
% integral loop will be costly. To improve performance, use parfor and
% compute pxpxi...for 4 gauss points in elements' loop, and take rows from
% them in inner loop.

% J=[px/pxi,py/pxi;px/peta,py/peta];
% px/pxi=sum(pN_i/pxi*x_i),etc.
x_i=elemNodePos(iElem,1:elemCfg.phyDOF:end)';
y_i=elemNodePos(iElem,2:elemCfg.phyDOF:end)';
if elemCfg.phyDOF==3
    z_i=elemNodePos(iElem,3:elemCfg.phyDOF:end)';
end
pxpxi=pNpxi*x_i;
pypxi=pNpxi*y_i;
pxpeta=pNpeta*x_i;
pypeta=pNpeta*y_i;
if elemCfg.phyDOF==3
    pzpxi=pNpxi*z_i;
    pzpeta=pNpeta*z_i;
    
    pxpzeta=pNpzeta*x_i;
    pypzeta=pNpzeta*y_i;
    pzpzeta=pNpzeta*z_i;
end
J_quad=zeros(elemCfg.phyDOF*elemCfg.phyDOF,elemCfg.nGaussPoint);
if elemCfg.phyDOF==2
    J_quad(1,:)=pxpxi(:);
    J_quad(2,:)=pxpeta(:);
    J_quad(3,:)=pypxi(:);
    J_quad(4,:)=pypeta(:);
elseif elemCfg.phyDOF==3
    J_quad(1,:)=pxpxi(:);
    J_quad(2,:)=pxpeta(:);
    J_quad(3,:)=pxpzeta(:);
    J_quad(4,:)=pypxi(:);
    J_quad(5,:)=pypeta(:);
    J_quad(6,:)=pypzeta(:);
    J_quad(7,:)=pzpxi(:);
    J_quad(8,:)=pzpeta(:);
    J_quad(9,:)=pzpzeta(:);
end