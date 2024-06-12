function [canRefine]=QC_can_be_refined(elemCfg,elemContainNetElem)
canRefine=zeros(elemCfg.nElem,1);
for iElem=1:elemCfg.nElem
    netElems=cell2mat(elemContainNetElem(:,1));
    if length(netElems)==3
        canRefine(iElem)=0;
    end
end
end