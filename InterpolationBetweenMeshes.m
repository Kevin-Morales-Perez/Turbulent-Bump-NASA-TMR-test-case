%Code to interpolate a solution from one grid to another with different
%size

%% 1.- Load solved grid variables 
%load("convergedRe93_61kturbulent.mat")%Turbulent flow at Re=93.61k
%load("converged_Re_2M_mesh8_turbulent.mat")%Turbulent flow at Re=1.0035 M
load("RampingSolmesh11_ItNum_35000.mat")
distPlat1=distPlat;


%% 2.- Add "_sol" prefix to u, v,p,nu_tilde fields and gradients to indicate
%are base fields

u_sol=u;
v_sol=v;
p_sol=p;
nu_tilde_sol=nu_tilde;

grad_u_sol=grad_u;
grad_v_sol=grad_v;
grad_p_sol=grad_p;
grad_nu_tilde_sol=grad_nu_tilde;

%% 3.-  Add sol suffix to cell center coordinates to indicate they correspond
%to the solved variable and store matrix sizes and number of elements

Xctrs_sol=Xctrs;
Yctrs_sol=Yctrs;
nx_sol=nx;
ny_sol=ny;
ncells_sol=ncells;

%% 4.- Generate new mesh
run("meshbump13_new_modified.m")

%% 5.- Process new mesh 

%New cells centroids
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

%Aditional Step
%This is to adjust the reference frames 
Xctrs_sol=Xctrs_sol + (distPlat - distPlat1);

%% 6.- Create a variable to store nearest point indexes from new mesh to
%solved mesh (indexes correspond to sol mesh)

indexes_new_to_sol=zeros(ny,nx,2);

%% 7.-Compute all the nearest points and store the indexes


%test
tic

for m=1:ny
    for n=1:nx


        x_i=Xctrs(m,n);
        y_i=Yctrs(m,n);

        %initialize minimum distance
        d_min=10000;



        for i=1:ny_sol
            for j=1:nx_sol
        
        
                %vector between central point and external points
                dist_vec=[Xctrs_sol(i,j),Yctrs_sol(i,j)]- [x_i,y_i];
        
                %get distance
                d_1=vecnorm(dist_vec);
                %compare
                d_min=min(d_1,d_min);
        
                %store indexes
                if d_min==d_1
        
                    %store coordinate
                    i_1=i;
                    j_1=j;
        
                    indexes_new_to_sol(m,n,1)=i_1;
                    indexes_new_to_sol(m,n,2)=j_1;
                end
        
            end
        end

    end
end

toc

%Test if are elements equal to zero in indexes_new_to_sol
for i=ny
    for j=nx

        a=indexes_new_to_sol(i,j,1);
        b=indexes_new_to_sol(i,j,2);

        if a==0 || b==0

            fprintf("Alert !")
            fprintf("j index = ")
            disp(a)
            fprintf("i index = ")
            disp(b)

        end 


        

    end
end 

%Success!!!


%% 8.-  Interpolate using gradients 
%Obtain u,v,p,nu tilde using an existing solution

%Init gradients 
grad_u_mn=zeros(1,2);
grad_v_mn=zeros(1,2);
grad_p_mn=zeros(1,2);
grad_nut_mn=zeros(1,2);

%INIT NEW FIELD VARIABLES
u=zeros(ny,nx); %velocity in X axis
v=zeros(ny,nx); %Velocity in Y axis
p=ones(ny,nx); %Pressure
nu_tilde=zeros(ny,nx);%EDDY Viscosity

for i=1:ny
    for j = 1:nx


        %Create vector from sol point to point in new grid 
        
        %new_point
        new_point = [Xctrs(i,j),Yctrs(i,j)];

        %Indexed of nearest point 
        m=indexes_new_to_sol(i,j,1);
        n=indexes_new_to_sol(i,j,2);

        %sol point 
        sol_point = [Xctrs_sol(m,n),Yctrs_sol(m,n)];

        vec_sol_new=new_point - sol_point;

        %U vel 
        grad_u_mn(1)=grad_u_sol(m,n,1);
        grad_u_mn(2)=grad_u_sol(m,n,2);
        
        u(i,j)=u_sol(m,n) + dot(grad_u_mn,vec_sol_new);

        %V vel
        grad_v_mn(1)=grad_v_sol(m,n,1);
        grad_v_mn(2)=grad_v_sol(m,n,2);
        
        v(i,j)=v_sol(m,n) + dot(grad_v_mn,vec_sol_new);


        %Press

        grad_p_mn(1)=grad_p_sol(m,n,1);
        grad_p_mn(2)=grad_p_sol(m,n,2);
        
        p(i,j)=p_sol(m,n) + dot(grad_p_mn,vec_sol_new);
        %Eddy Visc


        grad_nut_mn(1)=grad_nu_tilde_sol(m,n,1);
        grad_nut_mn(2)=grad_nu_tilde_sol(m,n,2);
        
        nu_tilde(i,j)=nu_tilde_sol(m,n) + dot(grad_nut_mn,vec_sol_new);

    end 
end 

%Success!!

save("NewInterpolationtomesh13.mat","u","v","p","nu_tilde")


%% Plot results

figure(4) %VELOCITY IN X AXIS
contourf(Xctrs,Yctrs,u, 20, 'LineColor', 'none')
title("Velocity in x axis (m/s)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

figure(5) %VELOCITY IN Y AXIS
contourf(Xctrs,Yctrs,v, 20, 'LineColor', 'none')
title("Velocity in y axis (m/s)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal


figure(8) %pressure
contourf(Xctrs,Yctrs,p, 20, 'LineColor', 'none')
title("Pressure (N/m2)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal


figure(13)%Nu tilde
contourf(Xctrs,Yctrs,nu_tilde, 20, 'LineColor', 'none')
title("Spalart - Allmaras nu~ transported variable")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal