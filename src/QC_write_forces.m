function force1=QC_write_forces(out1,boundaryCfg,materialCfg,projectCfg,result,iLoad,totalDisp)
    force = find_reaction_traction(boundaryCfg,result);
    kB=materialCfg.kB;
    T=materialCfg.T;
    b=materialCfg.b;
    if iLoad==1
        fprintf(out1,'iLoad lbd');
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
    fprintf(out1,'%d,%14.6e',iLoad,totalDisp/boundaryCfg.sz(2)+1);
    for i=1:length(force)
        fprintf(out1,',%14.6e',force(i)*b/kB/T);
    end
    if projectCfg.surfaceTraction
        for i=1:length(result.traction_forces)
            fprintf(out1,',%14.6e',result.traction_forces(i)*b/kB/T);
        end
    end
    fprintf(out1,'\n');
    force1=force(1)*b/kB/T;
end

