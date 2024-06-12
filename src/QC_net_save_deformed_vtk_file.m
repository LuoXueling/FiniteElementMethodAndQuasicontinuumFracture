function QC_net_save_deformed_vtk_file(datpath,iLoad,nodeCoordRef,elemNodeNo,elemCfg,projectCfg,disp_n1,velo_n1,acce_n1,result,QCbrokenLink)
% Write .vtk file for Paraview.
%Input the information
disp=zeros(length(nodeCoordRef),1);
try
%     disp(2:2:end)=result.disp_net(:,2)/result.totalDisp;
    disp(2:2:end)=result.disp_net(:,2);
    disp(1:2:end)=result.disp_net(:,1);
catch
    disp=result.disp_net;
%     disp(2:2:end)=disp(2:2:end)/result.totalDisp;
end
deformed_nodeCoordRef=nodeCoordRef+disp;
%Open file
fname=strcat(datpath,'/','QC_net_deformed_time_',num2str(iLoad),'.vtk');
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
iconst1=(elemCfg.nElem-length(find(QCbrokenLink)))*(elemCfg.nElemNode+1);

fprintf(out,'CELLS %5d  %5d  \n', elemCfg.nElem-length(find(QCbrokenLink)),iconst1);

for ielem=1:elemCfg.nElem
    if QCbrokenLink(ielem)==1
        continue;
    end
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

fprintf(out,'CELL_TYPES %5d\n',elemCfg.nElem-length(find(QCbrokenLink)));

for i=1:elemCfg.nElem
    if QCbrokenLink(i)==1
        continue;
    end
fprintf(out,'%2d\n',ntype);
end


end%if