function write_forces(out1,boundaryCfg,projectCfg,result,iLoad,totalDisp)
    force = find_reaction_traction(boundaryCfg,result);
    if iLoad==1
        fprintf(out1,'iLoad totalDisp');
        for i=1:length(force)
            fprintf(out1,strcat(' Reaction-',num2str(i)));
        end
        if projectCfg.surfaceTraction
            for i=1:length(result.traction_forces)
                fprintf(out1,strcat(' Traction-',num2str(i)));
            end
        end
        fprintf(out1,'\n');
    end
    fprintf(out1,'%d,%14.6e',iLoad,totalDisp);
    for i=1:length(force)
        fprintf(out1,',%14.6e',force(i));
    end
    if projectCfg.surfaceTraction
        for i=1:length(result.traction_forces)
            fprintf(out1,',%14.6e',result.traction_forces(i));
        end
    end
    fprintf(out1,'\n');
end

