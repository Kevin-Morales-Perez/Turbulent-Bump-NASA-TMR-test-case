%Mass Flux Test



%% Mesh test creation
nx=3;
ny=3;


X=zeros(ny+1,nx+1);
Y=zeros(ny+1,nx+1);

for i=1:ny+1
    X(i,:)=linspace(0,9,ny+1);
end

Y(ny+1,:)=[0,1,2,3];

for i=ny:-1:1
    Y(i,:)=Y(i+1,:)+3;
end

%Visualization of the grid.
figure(1);
% Plot vertical lines (constant j)
plot(X, Y, 'Color', 'b', 'LineWidth', 0.5);
hold on;
% Plot horizontal lines (constant i)
plot(X', Y', 'Color', 'b', 'LineWidth', 0.5);
axis equal;
xlabel('X');
ylabel('Y');
title('Mesh ');


%% Mesh Process

Xctrs=zeros(ny,nx);                 %Cell center X axis coordinate
Yctrs=zeros(ny,nx);                 %Cell center Y axis coordinate

cellVertxs=zeros(ny,nx,4,2);        %Coordinates of vertexes of each cell
cellCentrs=zeros(ny,nx,2);          %Centroids of each cell

lgtFaces=zeros(ny,nx,4);            %Face Areas
faceCentrs=zeros(ny,nx,4,2);        %Coordinates of face center by cell
%[fCw;fCn;fCe;fCs]

cellVols=zeros(ny,nx);              %Volumes of each cell

uVecNormFaces=zeros(ny,nx,4,2);     %Unitary vectors normal to faces

uVecNeighbNods=zeros(ny,nx,4,2);    %Unitary vectors from central node to 
% nb nodes

distNeighbNods=zeros(ny,nx,4);      %Distances between nodes per cell

uVecParlFaces=zeros(ny,nx,4,2);     %unitary vectors paralel to each face 
% of the cell (All from W to E and from S to N)

wlsqOperator=zeros(ny,nx,2,4);      %weighted least squares operators
%Example of use
%wlsop=reshape(wlsqOperator(i,j,:,:),[2,4])
%(wlsop*dif_vec)' ([2,4][4,1])'=[1,2]

VecCentNodFaceCents=zeros(ny,nx,4,2); %vectors from cell center to each ...
% face center

VecCentVertx=zeros(ny,nx,4,2);      %vectors from cell center to each ...
% vertex of the cell

VecCentNbNods=zeros(ny,nx,4,2);      %Vectors from cell center to NB cell
% centers

weightDistFactors=zeros(ny,nx,4);    %Weight distance factors (cell center)

distMinWall=zeros(ny,nx);            %Minimun distance from node to the 
% nearest wall
distPlat=0;

%Mesh process   _________________________________________________________ *
[cellVertxs,cellCentrs,lgtFaces,uVecParlFaces,faceCentrs,...
    cellVols,uVecNormFaces,distNeighbNods,uVecNeighbNods,...
    wlsqOperator,VecCentNodFaceCents,VecCentVertx,VecCentNbNods,...
    weightDistFactors,Xctrs,Yctrs,distMinWall] =...
    mesh_geometrical_process(X,Y,nx,ny,cellVertxs,cellCentrs,...
    lgtFaces,uVecParlFaces,faceCentrs,cellVols,uVecNormFaces,...
    distNeighbNods,uVecNeighbNods,wlsqOperator,VecCentNodFaceCents,...
    VecCentVertx,VecCentNbNods,weightDistFactors,Xctrs,Yctrs,distMinWall,...
    distPlat);



%% Field Variables

u=zeros(ny,nx); %velocity in X axis
v=zeros(ny,nx); %Velocity in Y axis
p=ones(ny,nx); %Pressure

p_prime=zeros(ny,nx);% Pressure correction

%Velocities normal to faces
u_f=zeros(ny,nx+1); %x velocity at faces
v_f=zeros(ny+1,nx); %y velocity at faces
