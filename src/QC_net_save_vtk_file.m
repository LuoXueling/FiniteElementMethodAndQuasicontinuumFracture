function QC_net_save_vtk_file(datpath,iLoad,nodeCoordRef,elemNodeNo,elemCfg,projectCfg,disp_n1,velo_n1,acce_n1,result,QCbrokenLink)
% Write .vtk file for Paraview.
%Input the information

deformed_nodeCoordRef=nodeCoordRef;
%Open file
fname=strcat(datpath,'/','QC_net_time_',num2str(iLoad),'.vtk');
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
        fprintf(out,'%14.6f  %14.6f  %14.6f \n',deformed_nodeCoordRef(2*ipoin-1),...
        deformed_nodeCoordRef(2*ipoin),0);
    end
elseif elemCfg.phyDOF==3
    for ipoin=1:elemCfg.nNode
        fprintf(out,'%14.6f  %14.6f  %14.6f \n',deformed_nodeCoordRef(3*ipoin-2),...
        deformed_nodeCoordRef(3*ipoin-1),deformed_nodeCoordRef(3*ipoin));
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
elseif(strcmp(elemCfg.elemType,'B21'))
    ntype=3;
end

fprintf(out,'CELL_TYPES %5d\n',elemCfg.nElem);

for i=1:elemCfg.nElem
fprintf(out,'%2d\n',ntype);
end

fprintf(out,'CELL_DATA %5d\n',elemCfg.nElem);
% 
%-------------------------------------------------------------------------
fprintf(out,'SCALARS fnet_normalized float  \n');
fprintf(out,'LOOKUP_TABLE default\n');

for ielem=1:elemCfg.nElem
    fprintf(out,'%14.6e\n',result.fnet(ielem));
end

fprintf(out,'SCALARS Broken_cell float  \n');
fprintf(out,'LOOKUP_TABLE default\n');
for iElem=1:elemCfg.nElem
    fprintf(out,'%d\n',QCbrokenLink(iElem));
end



end%if