function [netNodeInElem,netElemInElem,elemArea,elemContainNetNode,elemContainNetElem,canRefine]=QC_update_info(elemCfg,QCnodeCoord,QCelemNodeNo,nodeCoord,elemNodeNo,elemNodePos,projectCfg)

    %%%%%%%%%%%%%%%%%%%%%% Bound repnode 2 netnode %%%%%%%%%%%%%%%%%%%%%%%
    % This step should be applied after mesh refinement
%     elemNetNodeNo=zeros(elemCfg.nElem,elemCfg.nElemNode);
%     for iElem=1:elemCfg.nElem
%         for iNode=1:elemCfg.nElemNode
%             pos=elemNodePos(iElem,(iNode-1)*elemCfg.phyDOF+1:iNode*elemCfg.phyDOF);
%             [~,I]=pdist2(QCnodeCoord,pos,'euclidean','Smallest',1);
%             elemNetNodeNo(iElem,iNode)=I;
%         end
%     end
    %%%%%%%%%%%%%%%%%%%% End Bound repnode 2 netnode %%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%% Find netnode-link inside elem %%%%%%%%%%%%%%%%%%%%%
    % This step should be applied after mesh refinement
    netNodeInElem=zeros(length(QCnodeCoord),8);
    parfor iNode=1:length(QCnodeCoord)
        pos=QCnodeCoord(iNode,:);
        [~,I]=pdist2(nodeCoord,pos,'euclidean','Smallest',1);
        searchElem=find(sum(elemNodeNo==I,2));
        for iElem=1:8
            if iElem<=length(searchElem)
                xv=elemNodePos(searchElem(iElem),1:elemCfg.phyDOF:end);
                yv=elemNodePos(searchElem(iElem),2:elemCfg.phyDOF:end);
                xv=[xv,xv(1)];
                yv=[yv,yv(1)];
                isin=inpolygon(pos(1),pos(2),xv,yv);
                if isin
                    netNodeInElem(iNode,iElem)=searchElem(iElem);
                end
            else
                break;
            end
        end
    end
    netElemInElem=zeros(length(QCelemNodeNo),2);
    parfor iElem=1:length(QCelemNodeNo)
        nodeNo=QCelemNodeNo(iElem,:);
        nodeInElem=netNodeInElem(nodeNo,:);
        tab=tabulate(nodeInElem(nodeInElem~=0));
        inElem=find(tab(:,2)==2);
        for i=1:2
            if i<=length(inElem)
                netElemInElem(iElem,i)=inElem(i);
            else
                break;
            end
        end
    end
    %%%%%%%%%%%%%%%% End Find netnode-link inside elem %%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%% Element area %%%%%%%%%%%%%%%%%%%%%
    elemArea=zeros(elemCfg.nElem,1);
    elemContainNetElem={};
    elemContainNetNode={};
    canRefine=ones(elemCfg.nElem,1);
    parfor iElem=1:elemCfg.nElem
        nodePos=elemNodePos(iElem,:);
        A=nodePos(1:2);
        B=nodePos(3:4);
        C=nodePos(5:6);
        AB=B-A;
        AC=C-A;
        elemArea(iElem)=1/2*abs(det([AB;AC]));
        netElems=union(find(netElemInElem(:,1)==iElem),find(netElemInElem(:,2)==iElem));
        elemContainNetElem(iElem,1)={netElems};
        if length(netElems)<projectCfg.cantRefineBelow
            canRefine(iElem)=0;
        end
        netNodes=find(sum(netNodeInElem==iElem,2));
        elemContainNetNode(iElem,1)={netNodes};
    end

    %%%%%%%%%%%%%%%% End Element area %%%%%%%%%%%%%%%%%%%

end