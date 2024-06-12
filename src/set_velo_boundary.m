function velo_half=set_velo_boundary(velo_half,boundaryCfg,velo)
% Enforce velocity boundary.
s=size(boundaryCfg.dispConstraintArray);
for i=1:s(1)
    dispConstraintArray=cell2mat(boundaryCfg.dispConstraintArray(i,1));
    velo_half(dispConstraintArray)=boundaryCfg.dispConsDir(i)*velo;
end