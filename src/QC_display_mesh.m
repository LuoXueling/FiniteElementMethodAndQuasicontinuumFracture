function QC_display_mesh(plotname,nodeCoord,elemNodeNo,QCnodeCoord,QCelemNodeNo,figpath,focusElem)
if nargin<7
    focusElem.nodeCoord=[];
    focusElem.elemNodeNo=[];
end
s1=size(nodeCoord);
s2=size(elemNodeNo);
nNode=s1(1);
nElemNode=s2(2);
nElem=s2(1);

fig=figure('visible','off');

L=max(nodeCoord(:,1))-min(nodeCoord(:,1));
H=max(nodeCoord(:,2))-min(nodeCoord(:,2));
set(fig,'Position',[100,100,floor(2*L/H*360),floor(320)]);

subplot(1,2,1)
axis equal;
axis off;

VisualCurrent2D(nodeCoord,elemNodeNo)
if ~isempty(focusElem.elemNodeNo)
    VisualCurrent2D(focusElem.nodeCoord,focusElem.elemNodeNo,'r')
end


s1=size(QCnodeCoord);
s2=size(QCelemNodeNo);
nNode=s1(1);
nElemNode=s2(2);
nElem=s2(1);

subplot(1,2,2)
axis equal;
axis off;
VisualNet2D(QCnodeCoord,QCelemNodeNo)

print(fig, strcat(figpath,'/QC_',plotname),'-djpeg','-r300');

end

function VisualCurrent2D(LXY,LE,facecolor)
if nargin==2
    facecolor=[0.9,0.9,0.9];
    X=LXY(:,1);Y=LXY(:,2);
    patch(X(LE)',Y(LE)',ones(size(LE')),'EdgeColor',[0 0.4470 0.7410],'FaceColor',facecolor);
else
    X=LXY(:,1);Y=LXY(:,2);
    patch(X(LE)',Y(LE)',ones(size(LE')),'EdgeColor',[0 0.4470 0.7410],'FaceColor',facecolor,'FaceAlpha',0.3);
end
end

function VisualNet2D(LXY,LE)
X=LXY(:,1);Y=LXY(:,2);
patch(X(LE)',Y(LE)',ones(size(LE')),'EdgeColor',[0 0.4470 0.7410],'Marker','o','MarkerFaceColor','k','MarkerSize',1);
end