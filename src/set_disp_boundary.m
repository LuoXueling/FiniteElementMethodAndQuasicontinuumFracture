function [disp_n,saveConsValue]=set_disp_boundary(disp_n,boundaryCfg,computeCfg)
% Enforce displacement boundary and record values of constrained nodes.
s=size(boundaryCfg.dispConstraintArray);
for i=1:s(1)
    dispConstraintArray=cell2mat(boundaryCfg.dispConstraintArray(i,1));
    disp_n(dispConstraintArray)=boundaryCfg.dispConsDir(i)*computeCfg.totalDisp;
end
saveConsValue=disp_n(boundaryCfg.dispConstraint);