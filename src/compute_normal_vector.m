function normal_vector=compute_normal_vector(coord)
    s=size(coord);
    if s(1)==2 %line
        L1=[coord(1,:)-coord(2,:),0];
        L2=[0,0,1];
        vector=cross(L1,L2);
        % Normalization for direction and length.
        if vector(2) ~= 0
            normal_vector=sign(vector(2))*vector/norm(vector);
        else
            normal_vector=sign(vector(1))*vector/norm(vector);
        end
        normal_vector=normal_vector(1:2);
    elseif s(1)==3 || s(1)==4
        L1=coord(1,:)-coord(2,:);
        L2=coord(2,:)-coord(3,:);
        vector=cross(L1,L2);
        if vector(3) ~= 0
            normal_vector=sign(vector(3))*vector/norm(vector);
        elseif vector(2) ~= 0
            normal_vector=sign(vector(2))*vector/norm(vector);
        else
            normal_vector=sign(vector(1))*vector/norm(vector);
        end
    end
end