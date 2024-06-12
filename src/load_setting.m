function boundaryCfg = load_setting(boundaryCfg,meshInfo,elemCfg,projectCfg,elemNodeNo,elemNodePos)
surfaceTractionArray={};
if projectCfg.surfaceTraction
    DOF=elemCfg.phyDOF;
    disp('Using surface tension from .inp file.');
    Stsize=size(meshInfo.Stsets);
    StNo=Stsize(1);
    SfsetNames=meshInfo.Sfsets(:,1);
    for i=1:StNo
        if DOF==2
            tmp=cell2mat(meshInfo.Stsets(i,4));
            surfaceTractionArray(i,1)={tmp(1:2)};
        elseif DOF==3
            surfaceTractionArray(i,1)=meshInfo.Stsets(i,4);
        end
        SfID=find(strcmp(SfsetNames,cell2mat(meshInfo.Stsets(i,3))));
        surfaceTractionArray(i,2)=meshInfo.Sfsets(SfID,2);
    end
    % In abaqus, surface traction is defined on elements rather than faces
    % or edges. https://abaqus.uclouvain.be/English/SIMACAEMODRefMap/simamod-c-deformablesurf.htm
    % N'*P should be integrated on edges, so a scheme to find the nodes on
    % edges is necessary.
    s=size(surfaceTractionArray);
    nFace=size(elemCfg.faces);nFace=nFace(1);
    for i=1:s(1)
        elems=cell2mat(surfaceTractionArray(i,2));
        elemConnectNode=elemNodeNo(elems(:),:);
        avgNormalVector=zeros(1,DOF);
        exposedFaceCount=0;
        % For each element, find out elements that share elemCfg.nFaceNode
        % with 
        for iElem=1:length(elems)
            elemConnectElem={};
            for iNode=1:elemCfg.nElemNode
                n=elemConnectNode(iElem,iNode);
                f=find(elemNodeNo==n);
                nodeConnectElem=f-floor((f-1)/elemCfg.nElem)*elemCfg.nElem;
                elemConnectElem(iNode,1)={nodeConnectElem};
            end
            elemConnectElemArray=cell2mat(elemConnectElem(:,1));
            tab=tabulate(elemConnectElemArray);
            shareFaceElem=tab(find(tab(:,2)==elemCfg.nFaceNode),1); %#ok<FNDSB>
            deleteFace=[];
            for j=1:length(shareFaceElem)
                shareNode=intersect(elemNodeNo(shareFaceElem(j),:),elemConnectNode(iElem,:));
                for k=1:length(shareNode)
                    idx(k)=find(elemConnectNode(iElem,:)==shareNode(k));
                end
                for k=1:nFace
                    [A,~]=sort(elemCfg.faces(k,:));
                    [B,~]=sort(idx);
                    if sum(abs(A-B))==0
                        deleteFace=[deleteFace,k];
                        break;
                    end
                end
            end
            exposedFace=setdiff(1:nFace,deleteFace);
            exposedFaceNode=elemCfg.faces(exposedFace,:);
            
            nExposedFace=size(exposedFaceNode);
            
            for iFace=1:nExposedFace(1)
                nodes=exposedFaceNode(iFace,:);
                coord=zeros(length(nodes),DOF);
                for iNode=1:length(nodes)
                    node=nodes(iNode);
                    coord(iNode,:)=elemNodePos(elems(iElem),(node*DOF-(DOF-1)):(node*DOF));
                end
                avgNormalVector=avgNormalVector+compute_normal_vector(coord);
                exposedFaceCount=exposedFaceCount+1;
            end

            surfaceTractionArray(i,2+iElem)={exposedFaceNode};
        end
        % Average normal vector
        avgNormalVector=avgNormalVector/exposedFaceCount;
        avgNormalVector=avgNormalVector/norm(avgNormalVector);
        
        for iElem=1:length(elems)
            exposedFaceNode=cell2mat(surfaceTractionArray(i,2+iElem));
            nExposedFace=size(exposedFaceNode);
            angles=zeros(nExposedFace(1),1);
            for iFace=1:nExposedFace(1)
                nodes=exposedFaceNode(iFace,:);
                coord=zeros(length(nodes),DOF);
                for iNode=1:length(nodes)
                    node=nodes(iNode);
                    coord(iNode,:)=elemNodePos(elems(iElem),(node*DOF-(DOF-1)):(node*DOF));
                end
                normalVector=compute_normal_vector(coord);
                angles(iFace)=acosd(dot(normalVector,avgNormalVector));
            end
            exposedFace=find(angles<35);
            if length(exposedFace)~=1
                exposedFace=find(angles==min(angles));
                exposedFace=exposedFace(1);
                % warning(cell2mat(strcat('Multiple exposed faces detected on element ',{32},num2str(elems(iElem)),'.')));
            end
            surfaceTractionArray(i,2+iElem)={exposedFaceNode(exposedFace,:)};
        end
    end
end
boundaryCfg.surfaceTractionArray=surfaceTractionArray;

