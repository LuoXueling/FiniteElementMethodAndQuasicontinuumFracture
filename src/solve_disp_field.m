function [disp_n1,velo_n1,acce_n1,result]=solve_disp_field(method,projectCfg,boundaryCfg,computeCfg,materialCfg,...
                elemCfg,elemNodeNo,elemNodePos,elemDofLabel,disp_n,velo_n,acce_n,N,pNpxi,pNpeta,pNpzeta,M,assem,display)
% Solve displacement field through implicit, explicit or data-driven method.

if(method==0&&strcmp(projectCfg.constitutiveModel,'linear'))
    [disp_n1,velo_n1,acce_n1,result]=solve_linear_static_disp(projectCfg,boundaryCfg,computeCfg,materialCfg,...
        elemCfg,elemNodeNo,elemNodePos,elemDofLabel,disp_n,velo_n,acce_n,N,pNpxi,pNpeta,pNpzeta,M,assem,display);

elseif(method==1&&strcmp(projectCfg.constitutiveModel,'linear'))
    [disp_n1,velo_n1,acce_n1,result]=solve_linear_dynamic_disp(projectCfg,boundaryCfg,computeCfg,materialCfg,...
        elemCfg,elemNodeNo,elemNodePos,elemDofLabel,disp_n,velo_n,acce_n,N,pNpxi,pNpeta,pNpzeta,M,assem,display);

elseif(method==0&&strcmp(projectCfg.constitutiveModel,'neoHookean'))
    [disp_n1,velo_n1,acce_n1,result]=solve_neoHookean_static_disp(projectCfg,boundaryCfg,computeCfg,materialCfg,...
        elemCfg,elemNodeNo,elemNodePos,elemDofLabel,disp_n,velo_n,acce_n,N,pNpxi,pNpeta,pNpzeta,M,assem,display);

elseif(method==1&&strcmp(projectCfg.constitutiveModel,'neoHookean'))
    [disp_n1,velo_n1,acce_n1,result]=solve_neoHookean_dynamic_disp(projectCfg,boundaryCfg,computeCfg,materialCfg,...
        elemCfg,elemNodeNo,elemNodePos,elemDofLabel,disp_n,velo_n,acce_n,N,pNpxi,pNpeta,pNpzeta,M,assem,display);
end