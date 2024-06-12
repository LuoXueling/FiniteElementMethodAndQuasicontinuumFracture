function [nodeCoord,elemNodeNo,elemNodePos,elemCfg,boundaryCfg,nodeCoordRef,elemDofLabel,assem,newNodeCoord,QCbrokenElem]=QC_mesh_refine(Eref,QCnodeCoord,QCelemNodeNo,QCbrokenElem,nodeCoord,elemNodeNo,elemNodePos,elemCfg,boundaryCfg,projectCfg,meshInfo)
    fprintf('Refine mesh. %d nodes, %d elements -> ',length(nodeCoord),length(elemNodeNo));
    delElem=[];
    newNodeCoord=[];
    while ~isempty(Eref)
%         [~,~,~,~,~,canRefine]=QC_update_info(elemCfg,QCnodeCoord,QCelemNodeNo,nodeCoord,elemNodeNo,elemNodePos);
        elemNo=Eref(end);
        Eref(end)=[];
        if ~isempty(find(delElem==elemNo))% || ~canRefine(elemNo)
            continue
        end
        
        % Track longest edge
        [E0,E1,E2]=track_longest_edge(elemNo,nodeCoord,elemNodeNo,elemNodePos);
        if E0~=E1
            % If the E1 and E2 are away from E0, E0 is added to the queue
            % again.
            Eref=[Eref,E0];
        end
        % divide elements
        delE1=false;
        delE2=false;
%         if canRefine(E1) || (~isempty(E2) && canRefine(E2))
%             if canRefine(E1)
                [tmpNodeCoord1,tmpElemNodeNo1,tmpElemNodePos1]=divide_element(E1,elemNodePos,nodeCoord,elemNodeNo);
        
                elemNodeNo(length(elemNodeNo)+1:length(elemNodeNo)+2,:)=tmpElemNodeNo1;
                elemNodePos(length(elemNodePos)+1:length(elemNodePos)+2,:)=tmpElemNodePos1;
                QCbrokenElem(length(QCbrokenElem)+1:length(QCbrokenElem)+2)=0;
                delE1=true;
%             end
            if ~isempty(E2)
%                 if canRefine(E2)
                    [tmpNodeCoord1,tmpElemNodeNo2,tmpElemNodePos2]=divide_element(E2,elemNodePos,nodeCoord,elemNodeNo);
                    elemNodeNo(length(elemNodeNo)+1:length(elemNodeNo)+2,:)=tmpElemNodeNo2;
                    elemNodePos(length(elemNodePos)+1:length(elemNodePos)+2,:)=tmpElemNodePos2;
                    QCbrokenElem(length(QCbrokenElem)+1:length(QCbrokenElem)+2)=0;
                    delE2=true;
%                 end
            end
            nodeCoord(length(nodeCoord)+1,:)=tmpNodeCoord1;
            newNodeCoord=[newNodeCoord;tmpNodeCoord1];
%         end

        Eref(Eref==E1)=[];
        if ~isempty(E2)
            Eref(Eref==E2)=[];
        end

        for iref=1:length(Eref)
            if Eref(iref)>E1 && (~isempty(E2) && Eref(iref)>E2) && delE1 && delE2
                Eref(iref) = Eref(iref)-2;
            elseif (Eref(iref)>E1 && delE1) || (~isempty(E2) && Eref(iref)>E2 && delE2)
                Eref(iref) = Eref(iref)-1;
            end
        end

        delElem=[];
        if delE1
            delElem=[delElem,E1];
        end
        if delE2
            delElem=[delElem,E2];
        end
        elemNodeNo(delElem,:)=[];
        elemNodePos(delElem,:)=[];
        QCbrokenElem(delElem)=[];

        elemCfg.nElem=length(elemNodeNo);

        elemCfg.nNode=length(nodeCoord);
        elemCfg.meshDispDOF=elemCfg.nNode*elemCfg.phyDOF;
        elemCfg.meshPhaseDOF=elemCfg.nNode;

    end
    boundaryCfg=boundary_setting(boundaryCfg.constraint_direction,nodeCoord,elemCfg,projectCfg,meshInfo);
    [nodeCoordRef,elemDofLabel,~,assem]=support_matrices(nodeCoord,elemNodeNo,elemCfg);
    fprintf('%d nodes, %d elements.\n',length(nodeCoord),length(elemNodeNo));
end

function [E0,E1,E2]=track_longest_edge(elemNo,nodeCoord,elemNodeNo,elemNodePos,E1shareWith,E0)
    if nargin==4
        E1shareWith=find_share_edge_with(elemNo,elemNodeNo,elemNodePos);
        E0=elemNo;
        if isempty(E1shareWith)
            % the element E0 does not share hypotenuse with others as one
            % of their edges
            E1=elemNo;
            E2=[];
            return
        end
    end
    E2shareWith=find_share_edge_with(E1shareWith,elemNodeNo,elemNodePos);
    if elemNo==E2shareWith
        % the element elemNo (E1) shares the hypotenuse with E1shareWith
        E1=elemNo;
        E2=E1shareWith;
    elseif isempty(E2shareWith)
        % the element E1shareWith is on the edge and does not share hypotenuse with others as one
        % of their edges
        E0=elemNo;
        E1=E1shareWith;
        E2=[];
    else
        % track E2shareWith until its hypotenuse is the hypotenuse of E1shareWith, or is on the edge of the
        % geometry
        [E0,E1,E2]=track_longest_edge(E1shareWith,nodeCoord,elemNodeNo,elemNodePos,E2shareWith,elemNo);
    end
end

function [lens]=elem_edge_length(elemNo,elemNodePos)
    nodePos=elemNodePos(elemNo,:);
    A=nodePos(1:2);
    B=nodePos(3:4);
    C=nodePos(5:6);
    lens=zeros(3,1);
    lens(1)=norm(B-A);
    lens(2)=norm(C-B);
    lens(3)=norm(A-C);
end

function [shareEdgeWith]=find_share_edge_with(elemNo,elemNodeNo,elemNodePos)
% find the element whose edge is the hypotenuse of elemNo.
    lens=elem_edge_length(elemNo,elemNodePos); %lens(1):AB-C, lens(2):BC-A, lens(3):AC-B
    maxlen=find(lens==max(lens));
    maxlenEdge=elemNodeNo(elemNo,(maxlen==1)*[1,2]+(maxlen==2)*[2,3]+(maxlen==3)*[3,1]);
    shareNode1=elemNodeNo==maxlenEdge(1);
    shareNode2=elemNodeNo==maxlenEdge(2);
    shareNodeWith=find(sum(shareNode1+shareNode2,2)==2);
    shareEdgeWith=setdiff(shareNodeWith,elemNo);
end

function [tmpNodeCoord,tmpElemNodeNo,tmpElemNodePos]=divide_element(E1,elemNodePos,nodeCoord,elemNodeNo)
    lens=elem_edge_length(E1,elemNodePos); %lens(1):AB-C, lens(2):BC-A, lens(3):AC-B
    len1=find(lens==max(lens));
    maxlenEdge1=elemNodeNo(E1,(len1==1)*[1,2]+(len1==2)*[2,3]+(len1==3)*[3,1]);
    oppNodeNo=elemNodeNo(E1,(len1==1)*3+(len1==2)*1+(len1==3)*2);
    midpoint=mean(nodeCoord(maxlenEdge1,:),1);
    nNode=length(nodeCoord);
    tmpNodeCoord=midpoint;
    nodeCoord(nNode+1,:)=tmpNodeCoord;
    divE1NodeNo=[nNode+1,oppNodeNo,maxlenEdge1(1)]; %anti-clockwise
    if ~is_anti_clock_wise(nodeCoord,divE1NodeNo)
        divE1NodeNo=[maxlenEdge1(1),oppNodeNo,nNode+1];
        divE2NodeNo=[nNode+1,oppNodeNo,maxlenEdge1(2)];
    else
        divE2NodeNo=[maxlenEdge1(2),oppNodeNo,nNode+1];
    end
    tmpElemNodeNo(1,:)=divE1NodeNo;
    tmpElemNodeNo(2,:)=divE2NodeNo;
    tmp=nodeCoord(divE1NodeNo,:)';
    tmpElemNodePos(1,:)=tmp(:);
    tmp=nodeCoord(divE2NodeNo,:)';
    tmpElemNodePos(2,:)=tmp(:);
end

function [isAntiClockWise]=is_anti_clock_wise(nodeCoord,nodeNo)
    pos=nodeCoord(nodeNo,:);
    AB=pos(2,:)-pos(1,:);
    AC=pos(3,:)-pos(1,:);
    A=[AB;AC];
    if det(A)<0
        isAntiClockWise=false;
    else
        isAntiClockWise=true;
    end
end
