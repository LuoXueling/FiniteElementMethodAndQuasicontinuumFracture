 function Bu=compute_Bu(J,iGauss,elemCfg,pNpxi,pNpeta,pNpzeta)
% Compute gradient shape function when integrating for displacement field.
% pNpX=[pN/px;pN/py]=J^-1 [pN/pxi;pN/peta];
% N is a row matrix
% Bdisp(3,8) for 3 strain contributions and 8 nodal displacements.
    if elemCfg.phyDOF==2
        pNpX=J\[pNpxi(iGauss,:);pNpeta(iGauss,:)];
        Bu=zeros(3,elemCfg.elemDispDOF);
        for iNode=1:elemCfg.nElemNode
            tmp=[pNpX(1,iNode),0;...
                 0,pNpX(2,iNode);...
                 pNpX(2,iNode),pNpX(1,iNode)];
            Bu(:,iNode*2-1:iNode*2)=tmp;
        end
    elseif elemCfg.phyDOF==3
        pNpX=J\[pNpxi(iGauss,:);pNpeta(iGauss,:);pNpzeta(iGauss,:)];
        Bu=zeros(6,elemCfg.elemDispDOF);
        for iNode=1:elemCfg.nElemNode
            tmp=[pNpX(1,iNode),             0,      0;...
                 0            ,pNpX(2,iNode) ,      0;...
                 0            ,             0,pNpX(3,iNode);...
                 pNpX(2,iNode),pNpX(1,iNode) ,      0;...
                 0            ,pNpX(3,iNode) ,pNpX(2,iNode);...
                 pNpX(3,iNode),             0,pNpX(1,iNode)];
            Bu(:,iNode*3-2:iNode*3)=tmp;
        end
    end
    Bu=sparse(Bu);

end