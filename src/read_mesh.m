function [nodeCoord,elemNodeNo,meshInfo]=read_mesh(projectCfg)
% Read .inp file from ABAQUS.
disp(strcat('Reading .inp file: ',projectCfg.filename));
file=fopen(strcat('source/',projectCfg.filename,'.inp'),'r');
i=0; j=0;
nodeCoord=[]; elemNodeNo=[];
while ~feof(file)
    line=fgetl(file);
    if(strcmp(line,'*Node'))
        %Nodes
        line=fgetl(file);
        while ~contains(line,'*Element, type=')
            i=i+1;
            nodeCoord(i,:)=str2num(line);
            line=fgetl(file);
        end
        %Elements
        meshInfo.elemType=line(strfind(line,'type=')+5:end);
        line=fgetl(file);
        while ~contains(line,'*')
            j=j+1;
            elemNodeNo(j,:)=str2num(line);
            line=fgetl(file);
        end
        % User sections
        partNsets={};
        partElsets={};
        sections={};
        if contains(line,'*Nset') || contains(line,'*Elset')
            while ~contains(line,'*End Part') && ~contains(line,'Section')
                if contains(line,'*Nset')
                    s=size(partNsets);
                    h=s(1)+1;
                    isgenerate=contains(line,'generate');
                    partNsets(h,1)={line(strfind(line,'nset=')+5:end)};
                    partNsets(h,2)={[]};
                    line=fgetl(file);
                    while ~contains(line,'*Nset') && ~contains(line,'*Elset') &&...
                            ~contains(line,'Section') && ~contains(line,'*End Part')
                        c=cell2mat(partNsets(h,2));
                        l=str2num(line);
                        if isgenerate
                            l=l(1):l(3):l(2);
                        end
                        partNsets(h,2)={[c,l]};
                        line=fgetl(file);
                    end
                elseif contains(line,'*Elset')
                    s=size(partElsets);
                    h=s(1)+1;
                    isgenerate=contains(line,'generate');
                    if isgenerate
                        partElsets(h,1)={line(strfind(line,'elset=')+6:strfind(line,', gen')-1)};
                    else
                        partElsets(h,1)={line(strfind(line,'elset=')+6:end)};
                    end
                    
                    partElsets(h,2)={[]};
                    line=fgetl(file);
                    while ~contains(line,'*Nset') && ~contains(line,'*Elset') &&...
                            ~contains(line,'Section') && ~contains(line,'*End Part')
                        c=cell2mat(partElsets(h,2));
                        l=str2num(line);
                        if isgenerate
                            l=l(1):l(3):l(2);
                        end
                        partElsets(h,2)={[c,l]};
                        line=fgetl(file);
                    end
                end
            end
            if contains(line,'Section')
                while ~contains(line,'*End Part')
                    if contains(line,'Section')
                        line=fgetl(file);
                        s=size(sections);
                        h=s(1)+1;
                        sections(h,1)={line(strfind(line,'elset=')+6:strfind(line,', m')-1)};
                        sections(h,2)={line(strfind(line,'material=')+9:end)};
                    end
                    line=fgetl(file);
                end
            end
        end
        meshInfo.partNsets=partNsets;
        meshInfo.partElsets=partElsets;
        meshInfo.sections=sections;
        while ~contains(line,'*End Instance')
            line=fgetl(file);
        end
        line=fgetl(file);line=fgetl(file);
        % Read assembly
        Nsets={};
        Elsets={};
        Sfsets={};
        while ~contains(line,'*End Assembly')
            if contains(line,'*Nset')
                s=size(Nsets);
                h=s(1)+1;
                isgenerate=contains(line,'generate');
                Nsets(h,1)={line(strfind(line,'nset=')+5:strfind(line,', in')-1)};
                Nsets(h,2)={[]};
                line=fgetl(file);
                while ~contains(line,'*Nset') && ~contains(line,'*Elset') &&...
                        ~contains(line,'*Surface') && ~contains(line,'*End Assembly')
                    c=cell2mat(Nsets(h,2));
                    l=str2num(line);
                    if isgenerate
                        l=l(1):l(3):l(2);
                    end
                    Nsets(h,2)={[c,l]};
                    line=fgetl(file);
                end
            elseif contains(line,'*Elset')
                s=size(Elsets);
                h=s(1)+1;
                isgenerate=contains(line,'generate');
                Elsets(h,1)={line(strfind(line,'elset=')+6:strfind(line,', in')-1)};
                Elsets(h,2)={[]};
                line=fgetl(file);
                while ~contains(line,'*Nset') && ~contains(line,'*Elset') &&...
                        ~contains(line,'*Surface') && ~contains(line,'*End Assembly')
                    c=cell2mat(Elsets(h,2));
                    l=str2num(line);
                    if isgenerate
                        l=l(1):l(3):l(2);
                    end
                    Elsets(h,2)={[c,l]};
                    line=fgetl(file);
                end
            elseif contains(line,'*Surface')
                s=size(Sfsets);
                h=s(1)+1;
                Sfsets(h,1)={line(strfind(line,'name=')+5:end)};
                Sfsets(h,2)={[]};
                line=fgetl(file);
                while ~contains(line,'*')
                    Elname=strsplit(line,', ');
                    ElsetID=find(strcmp(Elsets(:,1),str2mat(Elname(1,1))));
                    c=cell2mat(Sfsets(h,2));
                    Sfsets(h,2)={[c,cell2mat(Elsets(ElsetID,2))]};
                    line=fgetl(file);
                end
            end
        end
        meshInfo.Nsets=Nsets;
        meshInfo.Elsets=Elsets;
        meshInfo.Sfsets=Sfsets;

        materials={};
        % Material definition
        if projectCfg.inpDefineMat
            try
                while ~contains(line,'MATERIALS')
                    line=fgetl(file);
                end
                line=fgetl(file);
                line=fgetl(file);
                while ~contains(line,'**')
                    if contains(line,'*Material')
                        s=size(materials);
                        h=s(1)+1;
                        materials(h,1)={line(strfind(line,'name=')+5:end)};
                        line=fgetl(file);
                        cnt=0;
                        while ~contains(line,'**') && ~contains(line,'*Material')
                            line=fgetl(file);
                            cnt=cnt+1;
                            materials(h,cnt+1)={str2num(line)};
                            line=fgetl(file);
                        end
                    end
                end
            catch
                error('No material found in .inp file.');
            end
        end
        meshInfo.materials=materials;
        % Boundary
        if strcmp(projectCfg.boundaryCfgFrom,'inpfile')
            try
                while ~contains(line,'** BOUNDARY CONDITIONS') && ~contains(line,'** LOADS')
                    line=fgetl(file);
                    if contains(line,'** LOADS')
                        warning('No boundary condition found in .inp file. Rigidbody displacement may cause bad convergence.');
                        break;
                    end
                end
                if contains(line,'** BOUNDARY CONDITIONS')
                    Bdsets={};
                    line=fgetl(file);line=fgetl(file);
                    while ~strcmp(line,'** ')
                        s=size(Bdsets);
                        h=s(1)+1;
                        Bdsets(h,1)={line(strfind(line,'Name: ')+6:strfind(line,' Type:')-1)};
                        Bdsets(h,2)={line(strfind(line,'Type: ')+6:end)};
                        line=fgetl(file);line=fgetl(file);
                        tmp=strsplit(line,', ');
                        Bdsets(h,3)=tmp(1,1);
                        s=size(tmp);
                        d=s(2)-2;
                        Bdsets(h,4:4+d)=tmp(1,2:end);
                        line=fgetl(file);
                    end
                    meshInfo.Bdsets=Bdsets;
                    if projectCfg.surfaceTraction
                        try
                            while ~contains(line,'** LOADS')
                                line=fgetl(file);
                            end
                        catch
                            error('No surface tension found in .inp file.');
                        end
                    end
                end
                if contains(line,'** LOADS')
                    Stsets={};
                    line=fgetl(file);line=fgetl(file);
                    while ~strcmp(line,'** ')
                        s=size(Stsets);
                        h=s(1)+1;
                        Stsets(h,1)={line(strfind(line,'Name: ')+6:strfind(line,'   Type: ')-1)};
                        Stsets(h,2)={line(strfind(line,'Type: ')+6:end)};
                        line=fgetl(file);line=fgetl(file);
                        tmp=strsplit(line,', ');
                        Stsets(h,3)=tmp(1,1);
                        dir=[str2num(cell2mat(tmp(1,4))),str2num(cell2mat(tmp(1,5))),str2num(cell2mat(tmp(1,6)))];
                        Stsets(h,4)={dir};
                        line=fgetl(file);
                    end
                    meshInfo.Stsets=Stsets;
                end
            catch
                if strcmp(projectCfg.boundaryCfgFrom,'inpfile')
                    warning('No boundary condition found in .inp file. Rigidbody displacement may cause bad convergence.');
                end
                if projectCfg.surfaceTraction
                    error('No surface tension found in .inp file.');
                end
            end
        end
    end
end
nodeCoord(:,1)=[];
elemNodeNo(:,1)=[];
save(strcat('source/',projectCfg.filename),'nodeCoord','elemNodeNo','meshInfo');
disp('Save mesh to .mat file.');
end