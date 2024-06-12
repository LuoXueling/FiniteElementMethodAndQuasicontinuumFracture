function [QCnodeCoord,QCelemNodeNo,QCmeshInfo]=QC_generate_network(path,Lx,Ly,s,crossPoint)
    if nargin==1
        Lx=64;
        Ly=64;
        s=1;
        crossPoint=true;
    end
    exceptPoly=[-1,96.1;224,96.1;255.9,128;224,159.9;-1,159.9];
    Nx=floor(Lx/s);
    Ny=floor(Ly/s);
    nNode=(Nx+1)*(Ny+1);
    QCnodeCoord=zeros(nNode,2);
    A=kron(0:Nx,ones(Ny+1,1))';
    B=kron(0:Ny,ones(Nx+1,1));
    QCnodeCoord(:,1)=A(:)*s;
    QCnodeCoord(:,2)=B(:)*s;
    QCelemNodeNo=[];

    h=waitbar(0,'Genrating elements');
    if crossPoint
        for i=0:Nx-1
            for j=0:Ny-1
                crossNodeNo=length(QCnodeCoord)+1;
                crossNodePos=(QCnodeCoord(1+i+(Nx+1)*j,:)+QCnodeCoord(Nx+3+i+(Nx+1)*j,:))/2;
                QCnodeCoord(crossNodeNo,:)=crossNodePos;
                if i==0 && j==0
                    QCelemNodeNo=[QCelemNodeNo;...
                        1+i+(Nx+1)*j,2+i+(Nx+1)*j;...
                        1+i+(Nx+1)*j,Nx+2+i+(Nx+1)*j;...
                        Nx+2+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        2+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        1+i+(Nx+1)*j,crossNodeNo;...
                        crossNodeNo,Nx+3+i+(Nx+1)*j;...
                        2+i+(Nx+1)*j,crossNodeNo;...
                        crossNodeNo,Nx+2+i+(Nx+1)*j];
                elseif i==0 && j~=0
                    QCelemNodeNo=[QCelemNodeNo;...
                        1+i+(Nx+1)*j,Nx+2+i+(Nx+1)*j;...
                        Nx+2+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        2+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        1+i+(Nx+1)*j,crossNodeNo;...
                        crossNodeNo,Nx+3+i+(Nx+1)*j;...
                        2+i+(Nx+1)*j,crossNodeNo;...
                        crossNodeNo,Nx+2+i+(Nx+1)*j];
                elseif i~=0 && j==0
                    QCelemNodeNo=[QCelemNodeNo;...
                        1+i+(Nx+1)*j,2+i+(Nx+1)*j;...
                        Nx+2+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        2+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        1+i+(Nx+1)*j,crossNodeNo;...
                        crossNodeNo,Nx+3+i+(Nx+1)*j;...
                        2+i+(Nx+1)*j,crossNodeNo;...
                        crossNodeNo,Nx+2+i+(Nx+1)*j];
                elseif i~=0 && j~=0
                    QCelemNodeNo=[QCelemNodeNo;...
                        Nx+2+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        2+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        1+i+(Nx+1)*j,crossNodeNo;...
                        crossNodeNo,Nx+3+i+(Nx+1)*j;...
                        2+i+(Nx+1)*j,crossNodeNo;...
                        crossNodeNo,Nx+2+i+(Nx+1)*j];
                end
                if mod(i*Nx+j,floor(Nx*Ny/100))==0
                    waitbar((i*Nx+j)/Nx/Ny);
                end
            end
        end
    else
        for i=0:Nx-1
            for j=0:Ny-1
                if i==0 && j==0
                    QCelemNodeNo=[QCelemNodeNo;...
                        1+i+(Nx+1)*j,2+i+(Nx+1)*j;...
                        1+i+(Nx+1)*j,Nx+2+i+(Nx+1)*j;...
                        Nx+2+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        2+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        1+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        2+i+(Nx+1)*j,Nx+2+i+(Nx+1)*j];
                elseif i==0 && j~=0
                    QCelemNodeNo=[QCelemNodeNo;...
                        1+i+(Nx+1)*j,Nx+2+i+(Nx+1)*j;...
                        Nx+2+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        2+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        1+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        2+i+(Nx+1)*j,Nx+2+i+(Nx+1)*j];
                elseif i~=0 && j==0
                    QCelemNodeNo=[QCelemNodeNo;...
                        1+i+(Nx+1)*j,2+i+(Nx+1)*j;...
                        Nx+2+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        2+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        1+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        2+i+(Nx+1)*j,Nx+2+i+(Nx+1)*j];
                elseif i~=0 && j~=0
                    QCelemNodeNo=[QCelemNodeNo;...
                        Nx+2+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        2+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        1+i+(Nx+1)*j,Nx+3+i+(Nx+1)*j;...
                        2+i+(Nx+1)*j,Nx+2+i+(Nx+1)*j];
                end
                if mod(i*Nx+j,floor(Nx*Ny/100))==0
                    waitbar((i*Nx+j)/Nx/Ny);
                end
            end
        end
    end
    close(h);
    h=waitbar(0,'Deleting nodes and elems.');

    delNodes=[];
    len1=length(QCnodeCoord);
    for iNode=1:len1
%     x=intersect(find(QCnodeCoord(:,1)>exceptArea(1,1)),find(QCnodeCoord(:,1)<exceptArea(1,2)));
%     y=intersect(find(QCnodeCoord(:,2)>exceptArea(2,1)),find(QCnodeCoord(:,2)<exceptArea(2,2)));
%     delNodes=intersect(x,y);
        if inpolygon(QCnodeCoord(iNode,1),QCnodeCoord(iNode,2),exceptPoly(:,1),exceptPoly(:,2))
            delNodes=[delNodes,iNode];
        end
    end
    QCnodeCoord(delNodes,:)=[];
    elemNodeNoDelMat=zeros(size(QCelemNodeNo));
    
    delElemMarker=zeros(length(QCelemNodeNo),1);

    len0=length(delNodes);
    for iNode=1:len0
        delNodeNo=delNodes(iNode);
        elemNodeNoDelMat=elemNodeNoDelMat+(QCelemNodeNo>delNodeNo);

        delElemMarker(QCelemNodeNo(:,1)==delNodes(iNode))=1;
        delElemMarker(QCelemNodeNo(:,2)==delNodes(iNode))=1;
        if mod(iNode,floor(len0/100))
            waitbar(iNode/len0);
        end
    end
    
    QCelemNodeNo=QCelemNodeNo-elemNodeNoDelMat;
    QCelemNodeNo(logical(delElemMarker),:)=[];
    close(h);

    QCmeshInfo.elemType='B21';
    QCmeshInfo.partNsets={};
    QCmeshInfo.partElsets={};
    QCmeshInfo.sections={};
    QCmeshInfo.Nsets={};
    QCmeshInfo.Elsets={};
    QCmeshInfo.Sfsets={};
    
    f=fopen(path,'w');
    fprintf(f,'*Heading\n** Job name: network Model name: network\n** Generated by: MATLAB script\n');
    fprintf(f,'**\n** PARTS\n**\n*Part, name=Part-1\n*Node\n');
    h=waitbar(0,'Writing nodes.');
    len1=length(QCnodeCoord);
    for iNode=1:len1
        fprintf(f,'\t%d,\t%.1f,\t%.1f\n',iNode,QCnodeCoord(iNode,1),QCnodeCoord(iNode,2));
        if mod(iNode,floor(len1/100))==0
            waitbar(iNode/len1);
        end
    end
    close(h);
    fprintf(f,'*Element, type=B21\n');
    h=waitbar(0,'Writing elems.');
    len2=length(QCelemNodeNo);
    for iElem=1:len2
        fprintf(f,'\t%d,\t%d,\t%d\n',iElem,QCelemNodeNo(iElem,1),QCelemNodeNo(iElem,2));
        if mod(iElem,floor(len2/100))==0
            waitbar(iElem/len2);
        end
    end
    fprintf(f,'*End Part\n**  \n**\n** ASSEMBLY\n**\n*Assembly, name=Assembly\n**  \n*Instance, name=Part-1-1, part=Part-1\n*End Instance\n**  \n*End Assembly');
    close(h);
end