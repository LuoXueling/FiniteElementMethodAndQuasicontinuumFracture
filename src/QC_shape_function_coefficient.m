function QCshapeFuncCoeff=QC_shape_function_coefficient(elemNodePos)
    x1=elemNodePos(:,1);
    y1=elemNodePos(:,2);
    x2=elemNodePos(:,3);
    y2=elemNodePos(:,4);
    x3=elemNodePos(:,5);
    y3=elemNodePos(:,6);
    a1=x2.*y3-x3.*y2;   b1=y2-y3;   c1=x3-x2;
    a2=x3.*y1-x1.*y3;   b2=y3-y1;   c2=x1-x3;
    a3=x1.*y2-x2.*y1;   b3=y1-y2;   c3=x2-x1;
    QCshapeFuncCoeff=[a1,a2,a3,b1,b2,b3,c1,c2,c3];
end