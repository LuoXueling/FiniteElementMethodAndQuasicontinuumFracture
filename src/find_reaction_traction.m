function [force] = find_reaction_traction(boundaryCfg,result)
% Find nodal force on constrained faces.
s=size(boundaryCfg.dispConstraintArray);
force=[];
for i=1:s(1)
    dispConstraintArray=cell2mat(boundaryCfg.dispConstraintArray(i,1));
    force(i)=sum(result.reaction(dispConstraintArray));
end
