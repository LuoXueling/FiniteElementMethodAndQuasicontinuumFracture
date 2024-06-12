function [projectCfg,jobpath,figpath,datpath,nodeCoord,elemNodeNo,meshInfo,...
          elemCfg,computeCfg,materialCfg,boundaryCfg,nodeCoordRef,...
          elemDofLabel,elemNodePos,assem,N,pNpxi,pNpeta,pNpzeta,M] = configs(configFrom)
% Set up everything essential for the computation
%% Display time
diary off;
time = string(datetime('now','Format','yyyy-MM-dd HH-mm-ss'));
disp(strcat('Start on',{32},time));
%% Project setting
if ~strcmp(configFrom,'script')
    load(configFrom,'projectCfg','nodeCoord','elemNodeNo','meshInfo',...
          'elemCfg','computeCfg','materialCfg','boundaryCfg','nodeCoordRef',...
          'elemDofLabel','elemNodePos','assem','N','pNpxi','pNpeta','pNpzeta','M');
    [jobpath,figpath,datpath]=create_dir(projectCfg.filename,time);
    save(strcat(jobpath,'/','configs'),'projectCfg','nodeCoord','elemNodeNo','meshInfo',...
          'elemCfg','computeCfg','materialCfg','boundaryCfg','nodeCoordRef',...
          'elemDofLabel','elemNodePos','assem','N','pNpxi','pNpeta','pNpzeta','M');
else
projectCfg.filename='Semisphere2D';
% 'inpfile' or 'script', only support Displacement constraints and surface traction,
% The constraint value only represent direction.
projectCfg.boundaryCfgFrom='inpfile'; 
% Load (Surface tension only) from .inp file.
projectCfg.surfaceTraction=false;
projectCfg.inpDefineMat=false;

projectCfg.noPhase=true;
projectCfg.QC=false;
projectCfg.networkFilename='QC-2-net';
projectCfg.QCfullyDiscrete=false;
projectCfg.QCrefineLeftRight=false;
projectCfg.QCrefineAroundNotch=false;
projectCfg.cantRefineBelow=10;
projectCfg.meshRefine=true;
projectCfg.fracture=false;
projectCfg.critlbd=0.5;
projectCfg.critElemFail=0.5;

projectCfg.constitutiveModel='neoHookean';
projectCfg.plane_stress_strain='stress'; %no effect for 3D cases;

projectCfg.Helem=2;
%% Create directories
[jobpath,figpath,datpath]=create_dir(projectCfg.filename,time);
%% Read mesh
% When processing a .inp file for the first time, the program runs 'catch',
% otherwise loads the corresponding .mat file.
if strcmp(projectCfg.boundaryCfgFrom,'inpfile')
    [nodeCoord,elemNodeNo,meshInfo]=read_mesh(projectCfg);
else
    try
        load(strcat('source/',projectCfg.filename,'.mat'));
        disp(cell2mat(strcat('Loading .mat file ',{32},projectCfg.filename)));
    catch
        [nodeCoord,elemNodeNo,meshInfo]=read_mesh(projectCfg);
    end 
end
%% Element definition
[elemCfg]=element_info(nodeCoord,elemNodeNo,meshInfo,projectCfg);
%% Compute setting
computeCfg.totalStep=2000;
% solve method: 'direct' for mldivide,
% 'auto' to automatically choose method
% 'distributed' for mldivide with distributed arrays
% 'cg' for Conjugate Gradient Method (not recommended until fully developed)
% 'pcg' for Preconditioned Conjugate Gradient Method (not recommended until fully developed)
% 'gpu' for gpuArray, but is not recommended unless the gpu has enough memory and is powerful enough

% 'direct' is recommended for small Kdisp
% 'distributed' is recommended for large Kdisp and cores>=4
computeCfg.solveDispMethod='direct'; 

computeCfg.disp=5;

computeCfg.staticVelo=-1; %mm/s
computeCfg.dynamicVelo=1;
computeCfg.staticPvelo=8; %MPa/s
computeCfg.dynamicPvelo=1000;

computeCfg.staticDt=0.01;
computeCfg.dynamicDt=1e-6;
computeCfg.critNormRes=1E-6;

computeCfg.detDisp=computeCfg.staticVelo*computeCfg.staticDt;
% vtkStep : write .vtk file per vtkStep step
computeCfg.vtkStep=10;
computeCfg.vtkDeformedStep=10;
% plotStep : plot .jpg file per plotStep step
computeCfg.plotStep=10000000;
% method : 0 for quasi-static, 1 for dynamic
computeCfg.method=0;

computeCfg=choose_method(computeCfg,projectCfg,elemCfg);
%% Material setting
% K : elastic stiffness matrix
% Gc : critical energy release rate
% rho : density

materialCfg.rho=1.e-6;
materialCfg.E=2285.71;
materialCfg.v=0.4286;
materialCfg.Gc=1.5;

if strcmp(projectCfg.constitutiveModel,'neoHookean')
    materialCfg.mu_pos=800;
    materialCfg.ratio=1; %mu_positive/mu_negative == K_positive/K_negative
    materialCfg.K_pos=2000;
    materialCfg.Gc=1.5;
    materialCfg.bulk=materialCfg.K_pos;
    materialCfg.shear=materialCfg.mu_pos;
end

if projectCfg.inpDefineMat
    materialCfg.materials={};
    materials=meshInfo.materials;
    s=size(materials);
    
    for iMat=1:s(1)
        tmpMat.name=cell2mat(materials(iMat,1));
        mat=cell2mat(materials(iMat,2));
        tmpMat.E=mat(1);
        tmpMat.v=mat(2);
        tmpMat.Gc=mat(3);
        if s(2)==3
            mat=cell2mat(materials(iMat,3));
            tmpMat.a=mat(1);
            tmpMat.C=mat(2);
        end
        tmpMat=material_setting(projectCfg,elemCfg,tmpMat);
        tmpMat.rho=materialCfg.rho;
        materialCfg.materials(iMat,1)={tmpMat};
    end

    materialCfg.elemMatType=zeros(elemCfg.nElem,1);
    sections=meshInfo.sections;
    elsets=meshInfo.partElsets;
    sSec=size(sections);
    sSet=size(elsets);

    for iSec=1:sSec(1)
        elname=cell2mat(sections(iSec,1));
        matname=cell2mat(sections(iSec,2));
        for iSet=1:sSet(1)
            name=cell2mat(elsets(iSet,1));
            if strcmp(name,elname)
                for iMat=1:s(1)
                    name=cell2mat(materials(iMat,1));
                    if strcmp(name,matname)
                        materialCfg.elemMatType(cell2mat(elsets(iSet,2)))=iMat;
                        break;
                    end
                end
                break;
            end
        end
    end
else
    tmpMat=material_setting(projectCfg,elemCfg,materialCfg);
    materialCfg=struct;materialCfg.materials={};
    materialCfg.materials(1,1)={tmpMat};
    materialCfg.elemMatType=ones(elemCfg.nElem,1);
end

materialCfg.kB=1.3806E-20;
materialCfg.T=200;
materialCfg.b=1E-17;
materialCfg.n=100;
materialCfg.Lc_L=13;
materialCfg.Lc=15;

%% Boundary setting
constraint_direction='Y';
boundaryCfg=boundary_setting(constraint_direction,nodeCoord,elemCfg,projectCfg,meshInfo);
%% Support matrices
[nodeCoordRef,elemDofLabel,elemNodePos,assem]=support_matrices(nodeCoord,elemNodeNo,elemCfg);
%% Shape function
[N,pNpxi,pNpeta,pNpzeta]=shape_function(meshInfo);
%% Compute mass
[M]=compute_mass(elemCfg,elemNodePos,N,pNpxi,pNpeta,pNpzeta,materialCfg,assem);
%% Load setting
boundaryCfg = load_setting(boundaryCfg,meshInfo,elemCfg,projectCfg,elemNodeNo,elemNodePos);
%% Save configurations
save(strcat(jobpath,'/','configs'),'projectCfg','nodeCoord','elemNodeNo','meshInfo',...
          'elemCfg','computeCfg','materialCfg','boundaryCfg','nodeCoordRef',...
          'elemDofLabel','elemNodePos','assem','N','pNpxi','pNpeta','pNpzeta','M');
end
