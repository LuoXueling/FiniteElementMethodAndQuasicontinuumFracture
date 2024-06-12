function [f_face_quad,result]=integrate_face_force(projectCfg,boundaryCfg,computeCfg,...
                elemCfg,elemNodePos)

elemDispDOF=elemCfg.elemDispDOF;
DOF=elemCfg.phyDOF;
nElem=elemCfg.nElem;
nFaceNode=elemCfg.nFaceNode;
W=elemCfg.gaussWeight;
P=computeCfg.P;
%%%%%%%%%%%%%%%%% Integrate face force %%%%%%%%%%%%%%%%%%%
f_face_quad=zeros(elemDispDOF,nElem);
result.traction_forces=0;
if projectCfg.surfaceTraction
    nForce=size(boundaryCfg.surfaceTractionArray);
    traction_forces=zeros(nForce(1),1);
    last_force=0;
    for iForce=1:nForce(1)
        dir=cell2mat(boundaryCfg.surfaceTractionArray(iForce,1))';
        elems=cell2mat(boundaryCfg.surfaceTractionArray(iForce,2));
        faces=boundaryCfg.surfaceTractionArray(iForce,3:end);
        for iElem=1:length(elems)
            f_face_e=zeros(elemDispDOF,1);
            nodes=cell2mat(faces(iElem));
            coord=zeros(length(nodes),DOF);
            for k=1:length(nodes)
                node=nodes(k);
                coord(k,:)=elemNodePos(elems(iElem),(node*DOF-(DOF-1)):(node*DOF));
            end
            area=compute_area(coord);
            for k=1:length(nodes)
                node=nodes(k);
                f_face_e((node*DOF-(DOF-1)):(node*DOF))=f_face_e((node*DOF-(DOF-1)):(node*DOF))+dir*P*area/nFaceNode;
            end
            f_face_quad(:,elems(iElem))=f_face_quad(:,elems(iElem))+f_face_e;
        end
        this_force=sum(sum(abs(f_face_quad)));
        traction_forces(iForce)=abs(this_force-last_force);
        last_force=this_force;
    end
    result.traction_forces=traction_forces;
end
%%%%%%%%%%%%%%% End Integrate face force %%%%%%%%%%%%%%%%%
end

function area=compute_area(coord)
    s=size(coord);
    if s(1)==2 %line
        area=norm(coord(1,:)-coord(2,:));
    elseif s(1)==3
        area=0.5*norm(cross(coord(2,:)-coord(1,:),coord(3,:)-coord(1,:)));
    elseif s(1)==4
        L1=coord(1,:)-coord(3,:);
        L2=coord(2,:)-coord(4,:);
        m=norm(L1);
        n=norm(L2);
        alpha=acos(dot(L1,L2)/(m*n));
        area=0.5*m*n*sin(alpha);
    end
end