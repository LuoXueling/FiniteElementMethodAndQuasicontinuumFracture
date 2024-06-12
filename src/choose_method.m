function computeCfg=choose_method(computeCfg,projectCfg,elemCfg)
% Automatically choose method for solving linear system.
    if  strcmp(computeCfg.solveDispMethod,'auto')
        if elemCfg.meshDispDOF>100000 || (strcmp(elemCfg.elemType,'C3D8')&&elemCfg.meshDispDOF>10000)
            computeCfg.solveDispMethod='distributed';
        else
            computeCfg.solveDispMethod='direct';
        end
        if strcmp(projectCfg.constitutiveModel,'linear')
            if elemCfg.meshDispDOF>70000
                computeCfg.solveDispMethod='distributed';
            end
        end
    end
end

