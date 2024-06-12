function acce_n1=set_acce_boundary(acce_n1,boundaryCfg)
% Enforce accelaration boundary.
acce_n1(boundaryCfg.dispConstraint(:))=0;