function boundaryCfg=boundary_setting(constraint_direction,nodeCoord,elemCfg,projectCfg,meshInfo)
% Find constrained coordinates based on loading direction.
boundaryCfg.constraint_direction=constraint_direction;
if strcmp(projectCfg.boundaryCfgFrom,'inpfile')
    DOF=elemCfg.phyDOF;
    try
        Bdsize=size(meshInfo.Bdsets);
%         disp('Using boundary condition from .inp file.');
    catch
        Bdsize=0;
    end
    BdNo=Bdsize(1);
    NsetNames=meshInfo.Nsets(:,1);
    boundaryCfg.dispConstraint=[];
    boundaryCfg.dispConstraintArray={};
    boundaryCfg.dispConsDir=[];
    for bd=1:BdNo
        NsetID=find(strcmp(NsetNames,cell2mat(meshInfo.Bdsets(bd,3))));
        dir=str2num(cell2mat(meshInfo.Bdsets(bd,4)));
        tmp=DOF*cell2mat(meshInfo.Nsets(NsetID,2))-(DOF-dir);
        boundaryCfg.dispConstraint=[boundaryCfg.dispConstraint,tmp];
        boundaryCfg.dispConstraintArray(bd,1)={tmp};
        try
            boundaryCfg.dispConsDir(bd)=sign(str2num(cell2mat(meshInfo.Bdsets(bd,6))));
        catch
            boundaryCfg.dispConsDir(bd)=0;
        end
    end
    boundaryCfg.dispConstraint=unique(boundaryCfg.dispConstraint)';
    boundaryCfg.dispFree=setdiff(1:elemCfg.meshDispDOF,boundaryCfg.dispConstraint);
    
    boundaryCfg.argXmin=find(nodeCoord(:,1)==min(nodeCoord(:,1)));
    boundaryCfg.argXmax=find(nodeCoord(:,1)==max(nodeCoord(:,1)));
    boundaryCfg.argYmin=find(nodeCoord(:,2)==min(nodeCoord(:,2)));
    boundaryCfg.argYmax=find(nodeCoord(:,2)==max(nodeCoord(:,2)));
    
    boundaryCfg.leftBottom=intersect(boundaryCfg.argXmin,boundaryCfg.argYmin);
    boundaryCfg.rightBottom=intersect(boundaryCfg.argXmax,boundaryCfg.argYmin);
    boundaryCfg.leftTop=intersect(boundaryCfg.argXmin,boundaryCfg.argYmax);
    boundaryCfg.rightTop=intersect(boundaryCfg.argXmax,boundaryCfg.argYmax);
    
    try
        boundaryCfg.corner=[boundaryCfg.leftBottom(1),boundaryCfg.rightBottom(1),boundaryCfg.leftTop(1),boundaryCfg.rightTop(1)];
        boundaryCfg.sz=[max(nodeCoord(:,1))-min(nodeCoord(:,1)),max(nodeCoord(:,2))-min(nodeCoord(:,2))];
    catch
    end
    
else
%     disp('Using boundary condition from script');
    %% Find the boundary of square geometry and the boundary for loading
    % dispConstraint : rows for different constraints, 3*-2 for x, 3*-1 for y, 3* for z.
    % dispFree : label of free displacements
    % dispConsDir : rows for different constraints, 1 for +, -1 for -;
    boundaryCfg.argXmin=find(nodeCoord(:,1)==min(nodeCoord(:,1)));
    boundaryCfg.argXmax=find(nodeCoord(:,1)==max(nodeCoord(:,1)));
    boundaryCfg.argYmin=find(nodeCoord(:,2)==min(nodeCoord(:,2)));
    boundaryCfg.argYmax=find(nodeCoord(:,2)==max(nodeCoord(:,2)));
    if elemCfg.phyDOF==3
        boundaryCfg.argZmin=find(nodeCoord(:,3)==min(nodeCoord(:,3)));
        boundaryCfg.argZmax=find(nodeCoord(:,3)==max(nodeCoord(:,3)));
    end
    
    boundaryCfg.leftBottom=intersect(boundaryCfg.argXmin,boundaryCfg.argYmin);
    boundaryCfg.rightBottom=intersect(boundaryCfg.argXmax,boundaryCfg.argYmin);
    boundaryCfg.leftTop=intersect(boundaryCfg.argXmin,boundaryCfg.argYmax);
    boundaryCfg.rightTop=intersect(boundaryCfg.argXmax,boundaryCfg.argYmax);
    
    try
        boundaryCfg.corner=[boundaryCfg.leftBottom(1),boundaryCfg.rightBottom(1),boundaryCfg.leftTop(1),boundaryCfg.rightTop(1)];
        boundaryCfg.sz=[max(nodeCoord(:,1))-min(nodeCoord(:,1)),max(nodeCoord(:,2))-min(nodeCoord(:,2))];
    catch
    end

    DOF=elemCfg.phyDOF;
    if strcmp(constraint_direction,'X')
        boundaryCfg.dispConstraintArray(1,1)={DOF*boundaryCfg.argXmax-(DOF-1)};
        boundaryCfg.dispConstraintArray(2,1)={DOF*boundaryCfg.argXmin-(DOF-1)};
    elseif strcmp(constraint_direction,'Y') 
        boundaryCfg.dispConstraintArray(1,1)={DOF*boundaryCfg.argYmax-(DOF-2)};
        boundaryCfg.dispConstraintArray(2,1)={DOF*boundaryCfg.argYmin-(DOF-2)};
        boundaryCfg.dispConstraintArray(3,1)={DOF*boundaryCfg.argYmax-(DOF-1)};
        boundaryCfg.dispConstraintArray(4,1)={DOF*boundaryCfg.argYmin-(DOF-1)};
    elseif strcmp(constraint_direction,'Z') 
        boundaryCfg.dispConstraintArray(1,1)={DOF*boundaryCfg.argZmax-(DOF-3)};
        boundaryCfg.dispConstraintArray(2,1)={DOF*boundaryCfg.argZmin-(DOF-3)};
    end
    tmp1=cell2mat(boundaryCfg.dispConstraintArray(1,1));
    tmp2=cell2mat(boundaryCfg.dispConstraintArray(2,1));
    tmp3=cell2mat(boundaryCfg.dispConstraintArray(3,1));
    tmp4=cell2mat(boundaryCfg.dispConstraintArray(4,1));
    boundaryCfg.dispConstraint=unique([tmp1(:);tmp2(:);tmp3(:);tmp4(:)]);

    boundaryCfg.dispFree=setdiff(1:elemCfg.meshDispDOF,boundaryCfg.dispConstraint);
    boundaryCfg.dispConsDir=[1;0;0;0];
end
