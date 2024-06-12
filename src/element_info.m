function [elemCfg]=element_info(nodeCoord,elemNodeNo,meshInfo,projectCfg)
% Define elements based on the type of element.

% nElem : num of elements
% nElemNode : num of nodes of each element
% nNode : num of nodes of the mesh
% nodeDOF : DOF of each node
% phyDOF : dimension of the structure
% nGaussPoint : num of Gauss points for integral
% elemDispDOF : DOF of displacement of each element
% elemPhaseDOF : DOF of phase of each element
% meshDispDOF : DOF of displacement of the mesh
% meshPhaseDOF : DOF of phase of the mesh
% Helem : size of elements in ABAQUS
% dt : discretization for loading

% 添加高阶单元，需要处理好形函数（每个单元哪几个形函数），以及节点插值时积分点分配的问题
[elemCfg.nElem,elemCfg.nElemNode]=size(elemNodeNo);
[elemCfg.nNode,elemCfg.nodeDOF]=size(nodeCoord);
elemCfg.phyDOF=elemCfg.nodeDOF;

elemCfg.elemDispDOF=elemCfg.nodeDOF*elemCfg.nElemNode;
elemCfg.elemPhaseDOF=elemCfg.nElemNode;

elemCfg.meshDispDOF=elemCfg.nNode*elemCfg.phyDOF;
elemCfg.meshPhaseDOF=elemCfg.nNode;

elemCfg.Helem=projectCfg.Helem;
if strcmp(meshInfo.elemType,'CPS4R') || strcmp(meshInfo.elemType,'CPS4')
    % The type should be CPS4 without reduced integration, but previously
    % generated .inp file used the type CPS4R.
    elemCfg.elemType='CPS4R';

    elemCfg.nGaussPoint=4;
    elemCfg.gaussWeight=[1,1,1,1];
    
    elemCfg.nFaceNode=2;
    elemCfg.faces=[1,2;...
                   2,3;...
                   3,4;...
                   1,4];
elseif strcmp(meshInfo.elemType,'CPS3')
%     warning('Using linear triangle elements which may cause irreliable result.');
    elemCfg.elemType='CPS3';
    elemCfg.nGaussPoint=1;
    elemCfg.gaussWeight=[1/2];
    elemCfg.nFaceNode=2;
    elemCfg.faces=[1,2;...
                   2,3;...
                   1,3];
%     elemCfg.nGaussPoint=3;
%     elemCfg.gaussWeight=[1/6,1/6,1/6];
elseif strcmp(meshInfo.elemType,'C3D4')
%     warning('Using linear tetrahedral elements which may cause irreliable result.');
    elemCfg.elemType='C3D4';
    elemCfg.nGaussPoint=1;
    elemCfg.gaussWeight=[1/6];
    elemCfg.nFaceNode=3;
    elemCfg.faces=[1,2,3;...
                   2,3,4;...
                   1,3,4;...
                   1,2,4];
%     elemCfg.nGaussPoint=4;
%     elemCfg.gaussWeight=[1/24,1/24,1/24,1/24];
elseif strcmp(meshInfo.elemType,'C3D8')
    % We do not use C3D8R in which 'R' means reduced integration with only
    % one gauss point that may cause mistakes in stress concentration.
    elemCfg.elemType='C3D8';
   
    elemCfg.nGaussPoint=8;
    elemCfg.gaussWeight=[1,1,1,1,1,1,1,1];
    elemCfg.nFaceNode=4;
    % clockwise or anticlockwise
    elemCfg.faces=[1,4,3,2;...
                   2,3,7,6;...
                   1,2,6,5;...
                   1,5,8,4;...
                   3,4,8,7;...
                   5,6,7,8]; 
elseif strcmp(meshInfo.elemType,'B21')
    elemCfg.elemType='B21';
   
    elemCfg.nGaussPoint=1;
    elemCfg.gaussWeight=[1];
    elemCfg.nFaceNode=1;

else
    error(strcat('Element not supported.',{32},meshInfo.elemType));
end
end