function [materialCfg]=material_setting(projectCfg,elemCfg,materialCfg)
% Set up material constants.
if strcmp(projectCfg.constitutiveModel,'linear')
    E=materialCfg.E; %Young's modulus
    v=materialCfg.v; %poisson's ratio
    K3D=E/(3*(1-2*v)); % ref: https://en.wikipedia.org/wiki/Bulk_modulus
    mu3D=E/(2*(1+v)); %shear modulus
    %ref: Eischen J W, Torquato S. Determining elastic behavior of
    %composites by the boundary element method[J]. Journal of applied physics, 1993, 74(1): 159-170.
    if strcmp(projectCfg.plane_stress_strain,'stress')
        K2D=9*K3D*mu3D/(3*K3D+4*mu3D);
        mu2D=mu3D;
    elseif strcmp(projectCfg.plane_stress_strain,'strain')
        K2D=K3D+mu3D/3;
        mu2D=mu3D;
    end
    % lame first param, ref:https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
    materialCfg.lambda=E*v/(1+v)/(1-2*v);
    % ref: https://en.wikipedia.org/wiki/Linear_elasticity
    if elemCfg.phyDOF==2
        materialCfg.KC=[K2D,K2D,0;
            K2D,K2D,0;
            0,  0,  0];
        materialCfg.muC=[mu2D,-mu2D,0;
            -mu2D,mu2D,0;
            0,   0    ,mu2D];
        materialCfg.K=materialCfg.muC+materialCfg.KC;
        materialCfg.shear=mu2D;
        materialCfg.bulk=K2D;
    elseif elemCfg.phyDOF==3
        materialCfg.KC=[K3D ,K3D ,K3D ,0   ,0  ,0;
            K3D ,K3D ,K3D ,0   ,0  ,0;
            K3D ,K3D ,K3D ,0   ,0  ,0;
            0   ,0   ,0   ,0   ,0  ,0;
            0   ,0   ,0   ,0   ,0  ,0;
            0   ,0   ,0   ,0   ,0  ,0];
        materialCfg.muC=[(4/3)*mu3D,-(2/3)*mu3D,-(2/3)*mu3D,0   ,0   ,0;
            -(2/3)*mu3D,(4/3)*mu3D,-(2/3)*mu3D,0   ,0   ,0;
            -(2/3)*mu3D,-(2/3)*mu3D,(4/3)*mu3D,0   ,0   ,0;
            0          ,0          ,0         ,mu3D,0   ,0;
            0          ,0          ,0         ,0   ,mu3D,0;
            0          ,0          ,0         ,0   ,0   ,mu3D];
        materialCfg.K=materialCfg.muC+materialCfg.KC;
        materialCfg.shear=mu3D;
        materialCfg.bulk=K3D;
    end
elseif strcmp(projectCfg.constitutiveModel,'neoHookean')
    materialCfg.K_neg=materialCfg.K_pos/materialCfg.ratio;
    materialCfg.mu_neg=materialCfg.mu_pos/materialCfg.ratio;
    if elemCfg.phyDOF==2
        K2D=materialCfg.K_pos;
        mu2D=materialCfg.mu_pos;
        materialCfg.KC=[K2D,K2D,0;
            K2D,K2D,0;
            0,  0,  0];
        materialCfg.muC=[mu2D,-mu2D,0;
            -mu2D,mu2D,0;
            0,   0    ,mu2D];
        materialCfg.K=materialCfg.muC+materialCfg.KC;
        materialCfg.shear=mu2D;
        materialCfg.bulk=K2D;
    elseif elemCfg.phyDOF==3
        K3D=materialCfg.K_pos;
        mu3D=materialCfg.mu_pos;
        materialCfg.KC=[K3D ,K3D ,K3D ,0   ,0  ,0;
            K3D ,K3D ,K3D ,0   ,0  ,0;
            K3D ,K3D ,K3D ,0   ,0  ,0;
            0   ,0   ,0   ,0   ,0  ,0;
            0   ,0   ,0   ,0   ,0  ,0;
            0   ,0   ,0   ,0   ,0  ,0];
        materialCfg.muC=[(4/3)*mu3D,-(2/3)*mu3D,-(2/3)*mu3D,0   ,0   ,0;
            -(2/3)*mu3D,(4/3)*mu3D,-(2/3)*mu3D,0   ,0   ,0;
            -(2/3)*mu3D,-(2/3)*mu3D,(4/3)*mu3D,0   ,0   ,0;
            0          ,0          ,0         ,mu3D,0   ,0;
            0          ,0          ,0         ,0   ,mu3D,0;
            0          ,0          ,0         ,0   ,0   ,mu3D];
        materialCfg.K=materialCfg.muC+materialCfg.KC;
        materialCfg.shear=mu3D;
        materialCfg.bulk=K3D;
    end
end
end