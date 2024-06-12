function Gu=compute_Gu(J,iGauss,elemCfg,pNpxi,pNpeta,pNpzeta)
% Compute displacement gradient shape function without Voigt transformation.
% pNpX=[pN/px;pN/py];

    if elemCfg.phyDOF==2
        pNpX=J\[pNpxi(iGauss,:);pNpeta(iGauss,:)];
        Gu=zeros(4,elemCfg.elemDispDOF);
        for iNode=1:elemCfg.nElemNode
            tmp=[pNpX(1,iNode),0;...
                 0,pNpX(2,iNode);...
                 pNpX(2,iNode),0;...
                 0,pNpX(1,iNode)];
            Gu(:,iNode*2-1:iNode*2)=tmp;
        end
    elseif elemCfg.phyDOF==3
        pNpX=J\[pNpxi(iGauss,:);pNpeta(iGauss,:);pNpzeta(iGauss,:)];
        Gu=zeros(9,elemCfg.elemDispDOF);
        for iNode=1:elemCfg.nElemNode
            tmp=[pNpX(1,iNode),             0,      0;...
                 0            ,pNpX(2,iNode) ,      0;...
                 0            ,             0,pNpX(3,iNode);...
                 pNpX(2,iNode),             0,      0;...
                 0            ,pNpX(3,iNode) ,      0;...
                 0            ,             0,pNpX(1,iNode);...
                 pNpX(3,iNode),             0,      0;...
                 0            ,pNpX(1,iNode) ,      0;...
                 0            ,             0,pNpX(2,iNode)];
            Gu(:,iNode*3-2:iNode*3)=tmp;
        end
    end
end