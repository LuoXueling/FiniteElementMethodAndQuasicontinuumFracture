function [nodeCoordRef,elemDofLabel,elemNodePos,assem]=support_matrices(nodeCoord,elemNodeNo,cfg)
% Compute matrices for convenience.
%% Vectorization
% nodeCoordRef : For 2D system, 2n for x, 2n+1 for y
nodeCoordRef=nodeCoord.';
nodeCoordRef=nodeCoordRef(:);
%% Find the label of DOF for each element
% Example : for 80000 elements, elemNodeNo is (80000x4), elemDofLabel is
% (80000x8) for 2D system.
elemDofLabel=zeros(cfg.nElem,cfg.nElemNode*cfg.nodeDOF);
for i=1:cfg.nodeDOF
    elemDofLabel(:,i:cfg.nodeDOF:end)=elemNodeNo*cfg.nodeDOF-(cfg.nodeDOF-i);
end
%% Find the position of nodes for each element
elemNodePos=nodeCoordRef(elemDofLabel);
sz=size(elemNodePos);
if sz(2)==1
    elemNodePos=elemNodePos';
end
%% Assemble matrices
% For a column in element, Me(4x1) for example, contributions will be
% distributed by elemNodeNo ((8x1) matrix by elemDofLabel).In parfor
% loops, assemble Me into Mquad(4xnElem), and use sparse(i,j,v,m,n) to
% assemble the system matrix (mxn). Function sparse can accumulate.
% Condition 1 : (4x1)
% Me(a)->M(elemNodeNo(iElem,a),1)
% Mquad(a,b)->M(elemNodeNo(b,a),1), so i=elemNodeNo.',j=1
assem.m4_1=elemNodeNo.';
assem.m4_1=assem.m4_1(:);
% Condition 2 : (8x1)
% Re(a)->R(elemDofLabel(iElem,a),1)
% Rquad(a,b)->R(elemDofLabel(b,a),1), so i=elemDofLabel.',j=1
assem.m8_1=elemDofLabel.';
assem.m8_1=assem.m8_1(:);
% Condition 3 : (4x4)
% Kphi_e(a,b)->Kphi(elemNodeNo(iElem,a),elemNodeNo(iElem,b))
% Kphi_quad(a,4*iElem+b)->Kphi(elemNodeNo(iElem,a),elemNodeNo(iElem,b))
% so i=kron(elemNodeNo.',ones(1,4)),j=kron(elemNodeNo.',ones(4,1))'
assem.m4_4i=kron(elemNodeNo.',ones(1,cfg.elemPhaseDOF));
assem.m4_4i=assem.m4_4i(:);
assem.m4_4j=kron(elemNodeNo.',ones(cfg.elemPhaseDOF,1));
assem.m4_4j=assem.m4_4j(:);
% Condition 4 : (8x8)
% Kdisp_e(a,b)->Kphi(elemDofLabel(iElem,a),elemDofLabel(iElem,b))
% Kdisp(a,8*iElem+b)->Kphi(elemDofLabel(iElem,a),elemDofLabel(iElem,b))
% so i=kron(elemNodeNo.',ones(1,8)),j=kron(elemNodeNo.',ones(8,1))';
assem.m8_8i=kron(elemDofLabel.',ones(1,cfg.elemDispDOF));
assem.m8_8i=assem.m8_8i(:);
assem.m8_8j=kron(elemDofLabel.',ones(cfg.elemDispDOF,1));
assem.m8_8j=assem.m8_8j(:);
% Condition 5 : (8x4)
assem.m8_4i=kron(elemDofLabel.',ones(1,cfg.elemPhaseDOF));
assem.m8_4j=kron(elemNodeNo.',ones(cfg.elemDispDOF,1));
% Condition 6 : (4x8)
assem.m4_8i=kron(elemNodeNo.',ones(1,cfg.elemDispDOF));
assem.m4_8j=kron(elemDofLabel.',ones(cfg.elemPhaseDOF,1));