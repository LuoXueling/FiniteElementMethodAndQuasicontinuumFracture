function solve(configFrom)
%% Configurations
[projectCfg,jobpath,figpath,datpath,nodeCoord,elemNodeNo,meshInfo,...
    elemCfg,computeCfg,materialCfg,boundaryCfg,nodeCoordRef,...
    elemDofLabel,elemNodePos,assem,N,pNpxi,pNpeta,pNpzeta,M] = configs(configFrom);
if projectCfg.QC
    QC_projectCfg = projectCfg;
    QC_projectCfg.filename = projectCfg.networkFilename;
    QC_projectCfg.boundaryCfgFrom = 'none';
    QC_projectCfg.inpDefineMat = false;
    QC_projectCfg.surfaceTraction = false;
    [QCnodeCoord,QCelemNodeNo,QCmeshInfo]=read_mesh(QC_projectCfg);
    [QCelemCfg]=element_info(QCnodeCoord,QCelemNodeNo,QCmeshInfo,QC_projectCfg);

    QC_display_mesh(num2str(0),nodeCoord,elemNodeNo,QCnodeCoord,QCelemNodeNo,figpath);

    if projectCfg.QCfullyDiscrete
        computeCfg.solveDispMethod='auto';
        computeCfg.solvePhaseMethod='auto';
    
        computeCfg=choose_method(computeCfg,QC_projectCfg,QCelemCfg);
    end
    [QCnodeCoordRef,QCelemDofLabel,QCelemNodePos,QCassem]=support_matrices(QCnodeCoord,QCelemNodeNo,QCelemCfg);
    QCboundaryCfg=boundary_setting('Y',QCnodeCoord,QCelemCfg,QC_projectCfg,QCmeshInfo);
    QCelemLength=zeros(length(QCelemNodeNo),1);
    parfor iElem=1:length(QCelemNodeNo)
        QCelemLength(iElem)=norm(QCnodeCoord(QCelemNodeNo(iElem,2),:)-QCnodeCoord(QCelemNodeNo(iElem,1),:));
    end
    QCbrokenLink=zeros(QCelemCfg.nElem,1);
    QCbrokenElem=zeros(elemCfg.nElem,1);

    QC_display_mesh(num2str(0),nodeCoord,elemNodeNo,QCnodeCoord,QCelemNodeNo,figpath);
    try
    if QC_projectCfg.QCrefineLeftRight
        for iRefine=1:7
            [netNodeInElem,~,~,~,~,~]=QC_update_info(elemCfg,QCnodeCoord,QCelemNodeNo,nodeCoord,elemNodeNo,elemNodePos,projectCfg);
            Eref=[];
            focusBoundary=[QCboundaryCfg.argXmin;QCboundaryCfg.argXmax];
            for iNetNode=1:length(focusBoundary)
                netNodeNo=focusBoundary(iNetNode);
                toberefine=netNodeInElem(netNodeNo,:);
                toberefine=toberefine(toberefine~=0);
                if length(toberefine)>1
                    continue;
                else
                    Eref=[Eref,toberefine];
                end
            end
            Eref=unique(Eref);
    %         Eref=1:elemCfg.nElem;
            [nodeCoord,elemNodeNo,elemNodePos,elemCfg,boundaryCfg,nodeCoordRef,elemDofLabel,assem,newNodeCoord,QCbrokenElem]=QC_mesh_refine(Eref,QCnodeCoord,QCelemNodeNo,QCbrokenElem,nodeCoord,elemNodeNo,elemNodePos,elemCfg,boundaryCfg,projectCfg,meshInfo);
            QC_display_mesh(num2str(iRefine),nodeCoord,elemNodeNo,QCnodeCoord,QCelemNodeNo,figpath);
        end
    end
    catch

    end

%     if projectCfg.QCfullyDiscrete
%         for iRefine=1:9
%             [netNodeInElem,~,~,~,~,~]=QC_update_info(elemCfg,QCnodeCoord,QCelemNodeNo,nodeCoord,elemNodeNo,elemNodePos,projectCfg);
%             Eref=1:1:elemCfg.nElem;
%             [nodeCoord,elemNodeNo,elemNodePos,elemCfg,boundaryCfg,nodeCoordRef,elemDofLabel,assem,newNodeCoord,QCbrokenElem]=QC_mesh_refine(Eref,QCnodeCoord,QCelemNodeNo,QCbrokenElem,nodeCoord,elemNodeNo,elemNodePos,elemCfg,boundaryCfg,projectCfg,meshInfo);
%             QC_display_mesh(num2str(iRefine),nodeCoord,elemNodeNo,QCnodeCoord,QCelemNodeNo,figpath);
%         end
%     end

    if QC_projectCfg.QCrefineAroundNotch
        [nodeCoord,elemNodeNo,elemNodePos,elemCfg,boundaryCfg,nodeCoordRef,elemDofLabel,assem,newNodeCoord,QCbrokenElem]=QC_mesh_refine(1:elemCfg.nElem,QCnodeCoord,QCelemNodeNo,QCbrokenElem,nodeCoord,elemNodeNo,elemNodePos,elemCfg,boundaryCfg,projectCfg,meshInfo);
        QC_display_mesh('0.0',nodeCoord,elemNodeNo,QCnodeCoord,QCelemNodeNo,figpath);
        for iRefine=1:3
            [netNodeInElem,~,~,~,~,canRefine]=QC_update_info(elemCfg,QCnodeCoord,QCelemNodeNo,nodeCoord,elemNodeNo,elemNodePos,projectCfg);
            Eref=[];
    
            aroundNotch=[intersect(find(QCnodeCoord(:,2)==160),find(QCnodeCoord(:,1)<224));...
                         intersect(find(QCnodeCoord(:,2)==96),find(QCnodeCoord(:,1)<224));...
                         intersect(intersect(find(QCnodeCoord(:,1)>=224),find(QCnodeCoord(:,1)<=256)),find(abs(QCnodeCoord(:,1)+QCnodeCoord(:,2)-384)<0.001));...
                         intersect(intersect(find(QCnodeCoord(:,1)>=224),find(QCnodeCoord(:,1)<=256)),find(abs(QCnodeCoord(:,1)-QCnodeCoord(:,2)-128)<0.001))];
    
            focusBoundary=[QCboundaryCfg.argXmin;QCboundaryCfg.argXmax;QCboundaryCfg.argYmin;QCboundaryCfg.argYmax;aroundNotch(:)];
            for iNetNode=1:length(focusBoundary)
                netNodeNo=focusBoundary(iNetNode);
                toberefine=netNodeInElem(netNodeNo,:);
                toberefine=toberefine(toberefine~=0);
                Eref=[Eref,toberefine];
            end
            Eref=unique(Eref);
    %         Eref=1:elemCfg.nElem;
            [nodeCoord,elemNodeNo,elemNodePos,elemCfg,boundaryCfg,nodeCoordRef,elemDofLabel,assem,newNodeCoord,QCbrokenElem]=QC_mesh_refine(Eref,QCnodeCoord,QCelemNodeNo,QCbrokenElem,nodeCoord,elemNodeNo,elemNodePos,elemCfg,boundaryCfg,projectCfg,meshInfo);
            QC_display_mesh(num2str(iRefine/10),nodeCoord,elemNodeNo,QCnodeCoord,QCelemNodeNo,figpath);
        end
    end
else
    QCbrokenLink=zeros(1,1);
    QCbrokenElem=zeros(1,1);
end

%% Set file object
outRT=fopen(strcat(datpath,'/','Reaction_Traction.dat'),'w');
out3=fopen(strcat(datpath,'/','driving_force.dat'),'w');
outener=fopen(strcat(datpath,'/','energy.dat'),'w');
%% Initialization
if projectCfg.QCfullyDiscrete
    [disp_n1,velo_n1,acce_n1,totalDisp,loadTime,P]=initialize(QCelemCfg);
else
    [disp_n1,velo_n1,acce_n1,totalDisp,loadTime,P]=initialize(elemCfg);
end

result=struct();
result.dDisp=sparse(elemCfg.meshDispDOF,1);
ext_work=0;
force=zeros(size(boundaryCfg.dispConstraintArray));
%
% totalTime=computeCfg.disp/computeCfg.staticVelo;
% syms iL coef
% timeFunc=computeCfg.staticDt*exp(-1/coef*iL);
% estimCoef=eval(vpasolve(int(timeFunc,iL,[1,computeCfg.totalStep])==totalTime,coef));
% fprintf('Estimate coeff = %.3f\n',estimCoef);
% timeFunc=subs(timeFunc,coef,estimCoef);
%% Main loop
for iLoad=1:computeCfg.totalStep
    tic
    %%%%%%%%%%%%%%%%%%%%%% Judge fracture condition %%%%%%%%%%%%%%%%%%%%%%
    %     dt=subs(timeFunc,iL,iLoad);
    % When crack initiate, start dynamic process
    if computeCfg.method==1 || ~isempty(find(QCbrokenLink))
        if projectCfg.QC
            method=0;
        else
            method=1;
           end
        dt=computeCfg.dynamicDt;
        Applied_velocity=computeCfg.dynamicVelo;
        totalDisp=totalDisp+Applied_velocity*dt;
        P=P+dt*computeCfg.dynamicPvelo;
        computeCfg.totalDisp=totalDisp;
    else
        method=0;
        dt=computeCfg.staticDt;
        Applied_velocity=computeCfg.staticVelo;
        totalDisp=totalDisp+Applied_velocity*dt;
        P=P+dt*computeCfg.staticPvelo;
        computeCfg.totalDisp=totalDisp;
    end
    %     if totalDisp-Applied_velocity*dt>computeCfg.disp
    %         break;
    %     end
    %%%%%%%%%%%%%%%%%%%% End Judge fracture condition %%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%% Update variables %%%%%%%%%%%%%%%%%%%%%%%%%%
    loadTime=loadTime+dt;
    computeCfg.detDisp=Applied_velocity*dt;
    computeCfg.P=P;
    computeCfg.iLoad=iLoad;
    fprintf('Step = %d/%d, time = %.3f, disp = %.3f\n',iLoad,computeCfg.totalStep,loadTime,totalDisp);
    %%%%%%% Updata field variables %%%%%%%%
    disp_n=disp_n1;velo_n=velo_n1;acce_n=acce_n1;
    %%%%%%%%%%%%%%%%%%%%%%%% End Update variables %%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Solve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%% Displacement field %%%%%%%%
    if projectCfg.QC
        if projectCfg.QCfullyDiscrete
            [disp_n1,result,QCbrokenLink]=QC_solve_discrete_disp(QC_projectCfg,QCboundaryCfg,computeCfg,materialCfg,...
                QCelemCfg,QCnodeCoord,QCelemNodeNo,QCelemNodePos,QCelemDofLabel,QCelemLength,disp_n,N,pNpxi,pNpeta,pNpzeta,M,QCassem,true,QCbrokenLink);
        else
            refineCnt=0;
            while true
                [disp_n1_tmp,disp_net,result,Eref,canRefine,QCbrokenLink,QCbrokenElem]=QC_solve_disp_field(projectCfg,boundaryCfg,computeCfg,materialCfg,...
                    elemCfg,nodeCoord,elemNodeNo,elemNodePos,elemDofLabel,QCelemLength,QCnodeCoord,QCelemNodeNo,QCelemDofLabel,QCelemCfg,QCboundaryCfg,QCbrokenLink,QCbrokenElem,disp_n,N,pNpxi,pNpeta,pNpzeta,M,assem,QCassem,true);
                if refineCnt>3
                    disp_n1=disp_n1_tmp;
                    break;
                end
                if ~isempty(Eref) && projectCfg.meshRefine
                    [nodeCoord,elemNodeNo,elemNodePos,elemCfg,boundaryCfg,nodeCoordRef,elemDofLabel,assem,newNodeCoord,QCbrokenElem]=QC_mesh_refine(Eref,QCnodeCoord,QCelemNodeNo,QCbrokenElem,nodeCoord,elemNodeNo,elemNodePos,elemCfg,boundaryCfg,projectCfg,meshInfo);
                    refineCnt=refineCnt+1;
                    if ~isempty(find(QCbrokenElem)) && projectCfg.fracture
                        focusElem.elemNodeNo=elemNodeNo(find(QCbrokenElem),:);
                        focusElem.nodeCoord=nodeCoord;
                        QC_display_mesh(strcat(num2str(iLoad),'_',num2str(refineCnt)),nodeCoord,elemNodeNo,QCnodeCoord,QCelemNodeNo,figpath,focusElem);
                    else
                        QC_display_mesh(strcat(num2str(iLoad),'_',num2str(refineCnt)),nodeCoord,elemNodeNo,QCnodeCoord,QCelemNodeNo,figpath);
                    end

                    s=size(newNodeCoord);
                    %%%%%%%% interpolate %%%%%%%
                    for iNewNode=1:s(1)
                        newNodePos=newNodeCoord(iNewNode,:);
                        [~,I]=pdist2(QCnodeCoord,newNodePos,'euclidean','Smallest',1);
                        disp_n1_tmp=[disp_n1_tmp;disp_net(I,1);disp_net(I,2)];
                    end
                    disp_n=disp_n1_tmp;
                    %%%%%%%% interpolate %%%%%%%
                else
                    disp_n1=disp_n1_tmp;
                    break;
                end

            end
        end
        velo_n1=velo_n;acce_n1=acce_n;
    else
        [disp_n1,velo_n1,acce_n1,result]=solve_disp_field(method,projectCfg,boundaryCfg,computeCfg,materialCfg,...
            elemCfg,elemNodeNo,elemNodePos,elemDofLabel,disp_n,velo_n,acce_n,N,pNpxi,pNpeta,pNpzeta,M,assem,true);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Solve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Optional actions %%%%%%%%%%%%%%%%%%%%%%%%%%
if(mod(iLoad,computeCfg.vtkStep)==0)
    disp('Save .vtk file');
    if projectCfg.QC
        QC_net_save_vtk_file(datpath,iLoad,QCnodeCoordRef,QCelemNodeNo,QCelemCfg,QC_projectCfg,disp_n1,velo_n1,acce_n1,result,QCbrokenLink);
        if ~projectCfg.QCfullyDiscrete
            QC_save_vtk_file(datpath,iLoad,nodeCoordRef,elemNodeNo,elemCfg,projectCfg,disp_n1,velo_n1,acce_n1,result,QCbrokenElem);
        end
    else
        save_vtk_file(datpath,iLoad,nodeCoordRef,elemNodeNo,elemCfg,projectCfg,disp_n1,velo_n1,acce_n1,result);
    end
end
if(mod(iLoad,computeCfg.vtkDeformedStep)==0)
    disp('Save deformed .vtk file');
    if projectCfg.QC
        result.totalDisp=totalDisp;
        QC_net_save_deformed_vtk_file(datpath,iLoad,QCnodeCoordRef,QCelemNodeNo,QCelemCfg,QC_projectCfg,disp_n1,velo_n1,acce_n1,result,QCbrokenLink);
        if ~projectCfg.QCfullyDiscrete
            QC_save_deformed_vtk_file(datpath,iLoad,nodeCoordRef,elemNodeNo,elemCfg,projectCfg,disp_n1,velo_n1,acce_n1,result,QCbrokenElem);
        end
    else
        save_deformed_vtk_file(datpath,iLoad,nodeCoordRef,elemNodeNo,elemCfg,projectCfg,disp_n1,velo_n1,acce_n1,result);
    end
end

% Reaction-Displacement
if projectCfg.QC
    if ~projectCfg.QCfullyDiscrete
        force=QC_write_forces(outRT,boundaryCfg,materialCfg,projectCfg,result,iLoad,totalDisp);
    else
        force=QC_write_forces(outRT,QCboundaryCfg,materialCfg,projectCfg,result,iLoad,totalDisp);
    end
    if force<1e-3 && iLoad>10
        disp('Force approx zero.');
        break;
    end
else
    write_forces(outRT,boundaryCfg,projectCfg,result,iLoad,totalDisp);
end

%%%%%%%%%%%%%%%%%%%%%%%% End Optional actions %%%%%%%%%%%%%%%%%%%%%%%%
toc
end
fclose('all');
diary off;
end

