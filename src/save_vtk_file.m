function save_vtk_file(datpath,iLoad,nodeCoordRef,elemNodeNo,elemCfg,projectCfg,disp_n1,velo_n1,acce_n1,result)
% Write .vtk file for Paraview.
%Input the information
sigmaXX=full(result.sigmaXX);
sigmaYY=full(result.sigmaYY);
sigmaXY=full(result.sigmaXY);

try
    SXX=full(result.SXX);
    SYY=full(result.SYY);
    SXY=full(result.SXY);
    PXX=full(result.PXX);
    PYY=full(result.PYY);
    PXY=full(result.PXY);
catch
end
epsilonXX=full(result.epsilonXX);
epsilonYY=full(result.epsilonYY);
epsilonXY=full(result.epsilonXY);
if elemCfg.phyDOF==3
    sigmaZZ=full(result.sigmaZZ);
    sigmaYZ=full(result.sigmaYZ);
    sigmaXZ=full(result.sigmaXZ);
    try
        SZZ=full(result.SZZ);
        SYZ=full(result.SYZ);
        SXZ=full(result.SXZ);
        PZZ=full(result.PZZ);
        PYZ=full(result.PYZ);
        PXZ=full(result.PXZ);
    catch
    end
    epsilonZZ=full(result.epsilonZZ);
    epsilonYZ=full(result.epsilonYZ);
    epsilonXZ=full(result.epsilonXZ);
end
disp=full(disp_n1);
velo=full(velo_n1);
acce=full(acce_n1);
%Open file
fname=strcat(datpath,'/','time_',num2str(iLoad),'.vtk');
out=fopen(fname,'w');

%Start writing:
%Header
fprintf(out,'# vtk DataFile Version 2.0\n');
fprintf(out,strcat('time_',num2str(iLoad),'.vtk\n'));
fprintf(out,'ASCII\n');
fprintf(out,'DATASET UNSTRUCTURED_GRID\n');

%Write nodal coordinates:
fprintf(out,'POINTS %5d float\n',elemCfg.nNode);

if elemCfg.phyDOF==2
    for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6f  %14.6f  %14.6f \n',nodeCoordRef(2*ipoin-1),...
        nodeCoordRef(2*ipoin),0);
    end
elseif elemCfg.phyDOF==3
    for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6f  %14.6f  %14.6f \n',nodeCoordRef(3*ipoin-2),...
        nodeCoordRef(3*ipoin-1),nodeCoordRef(3*ipoin));
    end
end

%Write element connectivity:
iconst1=elemCfg.nElem*(elemCfg.nElemNode+1);

fprintf(out,'CELLS %5d  %5d  \n', elemCfg.nElem,iconst1);

for ielem=1:elemCfg.nElem
    fprintf(out,'%5d  ',elemCfg.nElemNode);
    for inode=1:elemCfg.nElemNode
        fprintf(out,'%5d  ',(elemNodeNo(ielem,inode)-1));
    end
    fprintf(out,'\n');
end

%Write cell types:
% ref:https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
if(strcmp(elemCfg.elemType,'CPS4R'))
    ntype=9;
elseif(strcmp(elemCfg.elemType,'CPS3'))
    ntype=5;
elseif(strcmp(elemCfg.elemType,'C3D4'))
    ntype=10;
elseif(strcmp(elemCfg.elemType,'C3D8'))
    ntype=12;
end

fprintf(out,'CELL_TYPES %5d\n',elemCfg.nElem);

for i=1:elemCfg.nElem
fprintf(out,'%2d\n',ntype);
end

%Write Nodal scaler & vector values:
 fprintf(out,'POINT_DATA %5d\n',elemCfg.nNode);
% 

%-------------------------------------------------------------------------

if elemCfg.phyDOF==2
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS Displacement_X float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',disp(2*ipoin-1));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS Displacement_Y float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',disp(2*ipoin));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS Velocity_X float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',velo(2*ipoin-1));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS Velocity_Y float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',velo(2*ipoin));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS Accelation_X float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',acce(2*ipoin-1));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS Accelation_Y float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',acce(2*ipoin));
    end
    
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS sigmaXX float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',sigmaXX(ipoin));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS sigmaYY float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',sigmaYY(ipoin));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS sigmaXY float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',sigmaXY(ipoin));
    end
    try
        x=result.SXX;
        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS SXX float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',SXX(ipoin));
        end
        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS SYY float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',SYY(ipoin));
        end
        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS SXY float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',SXY(ipoin));
        end
        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS PXX float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',PXX(ipoin));
        end
        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS PYY float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',PYY(ipoin));
        end
        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS PXY float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',PXY(ipoin));
        end
    catch
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS epsilonXX float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',epsilonXX(ipoin));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS epsilonYY float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',epsilonYY(ipoin));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS epsilonXY float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',epsilonXY(ipoin));
    end
end
%-------------------------------------------------------------------------
if elemCfg.phyDOF==3
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS Displacement_X float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',disp(3*ipoin-2));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS Displacement_Y float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',disp(3*ipoin-1));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS Displacement_Z float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',disp(3*ipoin));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS Velocity_X float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',velo(3*ipoin-2));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS Velocity_Y float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',velo(3*ipoin-1));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS Velocity_Z float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',velo(3*ipoin));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS Accelation_X float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',acce(3*ipoin-2));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS Accelation_Y float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',acce(3*ipoin-1));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS Accelation_Z float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',acce(3*ipoin));
    end
    
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS sigmaXX float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',sigmaXX(ipoin));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS sigmaYY float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',sigmaYY(ipoin));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS sigmaZZ float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',sigmaZZ(ipoin));
    end
    
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS sigmaXY float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',sigmaXY(ipoin));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS sigmaYZ float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',sigmaYZ(ipoin));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS sigmaXZ float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',sigmaXZ(ipoin));
    end
    try
        x=result.SXX;
        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS SXX float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',SXX(ipoin));
        end
        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS SYY float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',SYY(ipoin));
        end
        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS SZZ float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',SZZ(ipoin));
        end

        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS SXY float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',SXY(ipoin));
        end
        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS SYZ float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',SYZ(ipoin));
        end
        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS SXZ float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',SXZ(ipoin));
        end
        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS PXX float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',PXX(ipoin));
        end
        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS PYY float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',PYY(ipoin));
        end
        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS PZZ float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',PZZ(ipoin));
        end

        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS PXY float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',PXY(ipoin));
        end
        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS PYZ float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',PYZ(ipoin));
        end
        %-------------------------------------------------------------------------
        fprintf(out,'SCALARS PXZ float  \n');
        fprintf(out,'LOOKUP_TABLE default\n');
        %
        for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6e\n',PXZ(ipoin));
        end
    catch
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS epsilonXX float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',epsilonXX(ipoin));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS epsilonYY float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',epsilonYY(ipoin));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS epsilonZZ float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',epsilonZZ(ipoin));
    end
    
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS epsilonXY float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',epsilonXY(ipoin));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS epsilonYZ float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',epsilonYZ(ipoin));
    end
    %-------------------------------------------------------------------------
    fprintf(out,'SCALARS epsilonXZ float  \n');
    fprintf(out,'LOOKUP_TABLE default\n');
    %
    for ipoin=1:elemCfg.nNode
    fprintf(out,'%14.6e\n',epsilonXZ(ipoin));
    end
end
fclose(out);
end%if