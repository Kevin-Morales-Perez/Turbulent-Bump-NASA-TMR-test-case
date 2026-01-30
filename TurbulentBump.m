%FINITE VOLUME METHOD FOR TURBULENT BUMP IN CHANNEL 
%TEST CASE FROM NASA TURBULENCE MODELLING RESOURCE
%https://turbmodels.larc.nasa.gov/bump.html
%ALGORITHM: SIMPLE
%CONVECTIVE SCHEME: UPWIND 1ST
% LINEAR SYST. SOL. METHOD: GAUSS-SEIDEL/GMRES-(ILU PRECONDITIONED)
%MASS FLUX INTERPOLATION METHOD: RHIE-CHOW
%CO-LOCATED NON-ORTHOGONAL STRUCTURED GRID
%STEADY STATE INCOMPRESSIBLE REYNOLDS AVERAGED NAVIER STOKES
%TURBULENCE MODEL: SPALART - ALLMARAS
%FLUID: AIR AT 20° AND SEA LEVEL 
%WRITTEN BY: KEVIN MORALES 
%AERONAUTICAL ENGINEERING - INSTITUTO POLITÉCNICO NACIONAL - CIUDAD DE
%MÉXICO

%CONTACT: 
%-kevin27182mora@gmail.com
%https://www.cfd-online.com/Forums/members/kevin+morales.html
%https://www.linkedin.com/in/kevin-morales-p%C3%A9rez-8928641a4/


%            BOUNDARY CONDITIONS FOR FLOW OVER A BUMP
%    ___________  NEUMMAN  (DU/DY=0 ,DV/DY=0)  _______________
%-->                                                       O N  --->
%--> I                                                     U E  ---> 
%--> N                                                     T U  DU/DX=0
%--> L                                                     L M  --->    
%--> E                                                     E M  DV/DX=0
%--> T                                                     T A  --->
%-->  ___ SYMETRIC ___|___   NO SLIP ___|___ SYMETRIC  ____  N  --->

close all
%clear all
%% DATA
load("D:\TurbulentBumpSecondVersionmfmod\RerampBoundTurbItNum_50000.mat")
%load("D:\TurbulentBumpSecondVersionmfmod\fourthSmoothingRe2M_mesh9_ItNum_5000.mat")
%load("D:\TurbulentBumpSecondVersionmfmod\finalConvergedRe3e6.mat")
%load("AcelConveg.mat")
%load precomputed fields
%load("interpolated_fields_Re300k_tomesh12.mat")
%load("finalConvergedRe3e6.mat")
%load("convergedRe3kmesh10_100fp.mat")
%load("RampingHighRe_mesh8_ItNum_50000.mat")
%load("converged_Re_1M_mesh8_turbulent.mat")
%load("interpolated_fields_Re1M.mat")%Reynolds 1 M
%load("convergedRe978_42kturbulent.mat")
%load("convergedRe93_61kturbulent.mat")
%load("interpolated_fields_Re93k.mat")
%physical parameters______________________________________________________*
rho=1.205;                    %Density (Kg/m3)
mu=1.802e-5;                  %Dynamic Viscosity (N*s/m^2)
nu = mu/rho;                  %Kinematic Viscosity (m/s^2)
%u0=0.045;                    %Velocity at the inlet (m/s) 
p0=1;                         %Outlet pressure (Prescribed)

%Geometrical parameters___________________________________________________*
bumpLgt=1.5;                %Length of the plate that contains the bump

%Reynolds Number _________________________________________________________*
Re=u0*rho*(bumpLgt-0.5)/mu; %Reynolds number

%Spalart Allmaras Model Constants_________________________________________*

%General constans
kappa=0.41; %Karman Constant
sigma_sa=2/3;%Turbulent Prantl Number used in SA model

%Basic constants
cb1=0.1355; 
cb2=0.622;

%Solid Damping
cv1=7.1;

%Wall Turbulence Destruction 
cw1=cb1/(kappa^2) +  (1 + cb2)/sigma_sa;
cw2=0.3;
cw3=2;

%% MESH 
%{
%Mesh generation _________________________________________________________*
%run MeshBump.m  to generate the mesh
%Geometry defined in MeshBump.m
%Comment once mesh generated if re-run
run("meshbump10_new_modified.m")
%run("meshbump6_new_modified.m")
%run("meshbump5_new_modified.m")
%run("meshbump4_new_modified.m")
%run("meshbump2_new_modified.m")
%run("meshbump2_new.m")
%run("meshbump2.m")

%Mesh variables _________________________________________________________*

%All element arrays are enumerated in each cell as follows
%[West,North,East,South], and Vertexces as [WN,EN,ES,WS]

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

%% FIELD VARIABLES

% Basic fields ___________________________________________________________

u=zeros(ny,nx); %velocity in X axis
v=zeros(ny,nx); %Velocity in Y axis
p=ones(ny,nx); %Pressure

p_prime=zeros(ny,nx);% Pressure correction

%Velocities normal to faces
u_f=zeros(ny,nx+1); %x velocity at faces
v_f=zeros(ny+1,nx); %y velocity at faces

u_f(:,1)=u0;%Imposing initial velocity at the inlet <----------------

%Vorticity
vorticity=zeros(ny,nx);
lower_omega=zeros(ny,nx);

%Gradient of velocity and pressure 
grad_u=zeros(ny,nx,1,2);
grad_v=zeros(ny,nx,1,2);
grad_p=zeros(ny,nx,1,2);
grad_p_prime=zeros(ny,nx,1,2);

%velocities at corners for cross - difussion  gradient computation
u_cor=zeros(ny,nx,1,4);%u
v_cor=zeros(ny,nx,1,4);%v

%   Turbulence fields_____________________________________________________

%Eddy Viscosity
mu_turbulent=zeros(ny,nx); %Eddy viscosity
nu_tilde=zeros(ny,nx) + 0.1*nu;%Kinematic Eddy modified viscosity + B.C

%Eddy viscosity at cell faces
mu_turbulent_fw=zeros(ny,nx);%Eddy viscosity at face W
mu_turbulent_fn=zeros(ny,nx);%Eddy viscosity at face N
mu_turbulent_fe=zeros(ny,nx);%Eddy viscosity at face E
mu_turbulent_fs=zeros(ny,nx);%Eddy viscosity at face S

%Nu~ at corners
nu_tilde_cor=zeros(ny,nx,1,4);

%Grad of Eddy Viscosity
grad_mu_turbulent=zeros(ny,nx,1,2);
%Grad of nu_tilde 
grad_nu_tilde=zeros(ny,nx,1,2);

% Reynolds stress tensor  Tau_ij
tau_xx=zeros(ny,nx); %-rho*u'2- Normal
tau_xy=zeros(ny,nx); %-rho*u'v'- Shear
tau_yy=zeros(ny,nx); %-rho*v'2- Normal

%Gradient of Reynolds Stresses
grad_tau_xx=zeros(ny,nx,1,2);
grad_tau_xy=zeros(ny,nx,1,2);
grad_tau_yy=zeros(ny,nx,1,2);

%%  COMPUTATION OF DIRECT AND NON ORTHOGONAL DIFFUSSION COEFFITIENTS
%Diffusive fluxes (Direct Gradient terms)
dW=zeros(ny,nx); %Face W
dN=zeros(ny,nx); %Face N
dE=zeros(ny,nx); %Face E
dS=zeros(ny,nx); %Face S

%Cross difussion terms (but they will be added as a source)
dW_c=zeros(ny,nx); %Face W
dN_c=zeros(ny,nx); %Face N
dE_c=zeros(ny,nx); %Face E
dS_c=zeros(ny,nx); %Face S

%compute coeffitients
[dW,dN,dE,dS,dW_c,dN_c,dE_c,dS_c] =diffusive_coeffitients(distNeighbNods,...
    uVecNormFaces,uVecNeighbNods,lgtFaces,uVecParlFaces,dW,dN,dE,dS,...
    dW_c,dN_c,dE_c,dS_c,nx,ny);
%% EQUATION SYSTEMS COEFFITIENTS 

%   Momentum equations ___________________________________________________
%Neighborhood coeffitients 
A_W=zeros(ny,nx);
A_N=zeros(ny,nx);
A_E=zeros(ny,nx);
A_S=zeros(ny,nx);

%Central coeffitient
A_P=zeros(ny,nx);
A_Pv=zeros(1,nx);%Diagonal Coeffitienst for simetric boundary condition
%at the free-stream zones

%Sources 
Su_x=zeros(ny,nx);
Su_y=zeros(ny,nx);

%   Pressure correction equations ________________________________________ 
%Neighborhood coeffitients 
Ap_W=zeros(ny,nx);
Ap_N=zeros(ny,nx);
Ap_E=zeros(ny,nx);
Ap_S=zeros(ny,nx);

%Central coeffitient
Ap_P=zeros(ny,nx);

%Sources 
Su_p=zeros(ny,nx);

%   Nu tilde equations ___________________________________________________
%Neighborhood coeffitients 
Ant_W=zeros(ny,nx);
Ant_N=zeros(ny,nx);
Ant_E=zeros(ny,nx);
Ant_S=zeros(ny,nx);

%Central coeffitient
Ant_P=zeros(ny,nx);

%Sources 
Su_nt=zeros(ny,nx);

%% INDEX CORRESPONDENCE BETWEEN MATRIX FORM AND VECTOR FORM

indx_mat=zeros(ncells,2);% MAtrix indexes
indx_k=zeros(ncells,1);% Vector indexes 

%k=(i-1)*nx + j

%P(center)  (i,j)  --->  k
%West       (i,j-1) ---> k - 1
%North      (i-1,j) ---> k - nx
%East       (i,j+1) ---> k + 1 
%South      (i+1,j) ---> k + 1

for i=1:ny
    for j=1:nx
        k=(i-1)*nx + j;
        indx_k(k)=k;
        indx_mat(k,:)=[i,j];
    end
end

%Interior, Edges , Corner Indexes

%vector indexes
indx_cor_WN=1;%West North corner (1)
indx_edg_N=2:nx-1;%North edge (2)
indx_cor_EN=nx;%East North corner(3)
indx_edg_W=nx+1:nx:(ny-2)*nx+1;%West Edge (4)
indx_edg_E=nx*2:nx:(ny-1)*nx;%East Edge (5)
indx_cor_WS=(ny-1)*nx +1;%West South corner (6)
indx_edg_S=(ny-1)*nx +2:ncells-1;%South Edge (7)
indx_cor_ES=ncells;%East South corner (8)
indx_interior=zeros(ny-2,nx-2);%indexes for interior cells

for i=2:ny-1
    indx_interior(i-1,:)=(i-1)*nx+2:1:i*nx-1;
end
indx_interior=reshape(transpose(indx_interior),1,[]);
%}
%% GAUSS SEIDEL SOLVER SETTINGS

%Underelaxation factors (Works also for GRMES)
alpha_uv=0.7;%0.1925;  %<---------------------     x-y momentum 
alpha_p=0.3; %<---------------------     pressure correction
alpha_nt=1; %<---------------------     nu-tilde

%Tolerance for inner iterations 
epsilon_u=1e-18;
epsilon_v=1e-18;
epsilon_p=1e-18;
epsilon_nt=1e-18;

%Max inner iterations
max_iterations_u=10;% Max iterations for momentum eqs.
max_iterations_v=10;% Max iterations for momentum  y eq .
max_iterations_p=50;% Max iterations for pressure eq.
max_iterations_nt=2;% Max iterations for SA equation.

%% GMRES/BICGSTAB SETTINGS
%x momentum
A_xsp=sparse(ncells,ncells);%ALL coeffitients
suX_vec=zeros(ncells,1);%Source Vector
restart_x=30;%Restart vector space
tol_x=1e-9;%tolerance desired
maxit_x=4;%Inner it

err_x2=1;%Output error from GRMES

%y momentum
A_ysp=sparse(ncells,ncells);%ALL coeffitients
suY_vec=zeros(ncells,1);%Source Vector
restart_y=30;%Restart vector space
tol_y=1e-9;%tolerance desired
maxit_y=4;%Inner it

err_y2=1;%Output error from GRMES

%pressure correction
A_psp=sparse(ncells,ncells);%ALL coeffitients
suP_vec=zeros(ncells,1);%Source Vector
restart_p=30;%Restart vector space
tol_p=1e-2;%tolerance desired
maxit_p=4;%Inner it

err_p2=1;%Output error from GRMES

%Spalart - Allmaras model equation
A_ntsp=sparse(ncells,ncells);%ALL coeffitients
suNt_vec=zeros(ncells,1);%Source Vector
restart_nt=2;%Restart vector space
tol_nt=1e-15;%tolerance desired
maxit_nt=3;%Inner it

err_nt2=1;%Output error from GRMES


%% MAIN SOLVER SETTINGS AND VARIABLES

%Solver for linear systems
%G-S:GAUSS SEIDEL
%GMRES:GENERALIZED MINIMUM RESIDUAL METHOD
%BICGSTAB:BICONJUGATE GRADIENT STABILIZED METHOD
solver_x_momentum="G-S";
solver_y_momentum="G-S";
solver_p_correction="G-S";
solver_nu_tilde="G-S";

%TVD scheme (DISABLED - ENABLED) 
tvd_mode="DISABLED";

%Impose again inlet conditions for precomputed fields
%u_f(:,1)=u0;%Imposing initial velocity at the inlet <----------------

max_iterations=300000;% Max outer iterations * <---------------------

error_tgt=1e-18; %Target error*
max_residual=1e10; % unstable value flag*
convergedFlg=false; %Flag for convergence*

iterations_trig=0;%Trigger conditions at certain iterations
ramping_cont=0;%Counter for ramping block
f_ramp=1;%.00011513588;%.00032899483;%Ramping factor
ramping_steps=0;% Number of steps to ramp
ramping_spacing=1;

%Residuals

%Raw Residuals
rsid_x=zeros(ny,nx); %X momentum*
rsid_y=zeros(ny,nx); %Y momentum*
rsid_p=zeros(ny,nx);% Pressure correction eq.*
rsid_nt=zeros(ny,nx);% Turbulence*
rsid_cont=zeros(ny,nx); % Continuity*

%Error from residual (L2 norms except continuity)*
err_x=1;%x momentum*
err_y=1;%y momentum*
err_p=1;%pressure correction*
err_cont=1;%continuity*
err_nt=1;%nu tilde*

residualsMat=zeros(max_iterations,5);% Residual Matrix*
inneritcontMat=zeros(max_iterations,4);% Inner iterations counter matrix* 
gmrsErrors=zeros(max_iterations,4);%Errors from GMRES if used*
turbulenceMonitoring=zeros(max_iterations,3);%Monitoring turbulence

%Inner iterations counters*
it_innx=0;
it_inny=0;
int_innp=0;
int_innnt=0;


%Outer iterations counter*
iterations_cont=0;

%% SOLVER

tic
%Main solver loop (Semi-Implicit algorithm for Pressure Linked Equations) 
while convergedFlg==false

    %Iterations Counter*
    iterations_cont=iterations_cont+1;

    %----------------------   I LAMINAR   --------------------------------

    % 0.1- compute pressure gradient *
    grad_p=computePressGradient(p,wlsqOperator,p0,f1p,f2p,nx,ny);

    % 0.2- compute velocity at corners for non orthogonal difussion *
    u_cor = computePhiVertex(u,grad_u,VecCentVertx,nx,ny);
    v_cor= computePhiVertex(v,grad_v,VecCentVertx,nx,ny);

    % 0.3.- compute eddy viscosity at cell faces *
    [mu_turbulent_fw,mu_turbulent_fn,mu_turbulent_fe,mu_turbulent_fs]...
        = computePhiFaces(mu_turbulent,weightDistFactors,nx,ny);

    switch tvd_mode

        case "DISABLED"
            % #1.- MOMENTUM LINK COEFFITIENTS
            [A_W,A_E,A_N,A_S,A_P,A_Pv,Su_x,Su_y] = momentum_link_coeff(rho,mu,...
                mu_turbulent_fw,mu_turbulent_fn,mu_turbulent_fe,mu_turbulent_fs...
                ,lgtFaces,cellVols,dW,dN,dE,dS,dW_c,dN_c,dE_c,dS_c,u0,u_f,v_f,...
                u_cor,v_cor,grad_p,nx,ny,solidMask);
        case "ENABLED"
            [aW,aE,aN,aS,aP,aPv,suX,suY] = momentum_link_coeffTVD(rho,mu,...
                mu_turbulent_fw,mu_turbulent_fn,mu_turbulent_fe,mu_turbulent_fs,...
                lgtFaces,cellVols,dW,dN,dE,dS,dW_c,dN_c,dE_c,dS_c,u0,u_f,v_f,...
                u_cor,v_cor,grad_p,nx,ny,solidMask,u,v,grad_u,grad_v,VecCentNbNods);
    end

    % #2.- SOLVE X MOMENTUM
    
    switch solver_x_momentum
        case "G-S"
    
            [u,rsid_x,err_x,it_innx] = GaussSeidel_GeneralSolver(A_P,A_W,...
            A_N,A_E,A_S,Su_x,rsid_x,epsilon_u,err_x,max_iterations_u,...
            alpha_uv,nx,ny,u,u);
        case "GMRES"
            [u,rsid_x,err_x,err_x2,it_innx] = gmresSolver(A_P,A_W,A_N,...
                A_E,A_S,Su_x,A_xsp,suX_vec,restart_x,tol_x,maxit_x,...
                indx_interior,indx_edg_W,indx_edg_N,indx_edg_E,indx_edg_S,...
                indx_cor_WN,indx_cor_EN,indx_cor_ES,indx_cor_WS, ...
                indx_mat,ny,nx,ncells,u,alpha_uv);

        case "BICGSTAB"
            [u,rsid_x,err_x,err_x2,it_innx] = bicgstabSolver(A_P,A_W,A_N,...
                A_E,A_S,Su_x,A_xsp,suX_vec,restart_x,tol_x,maxit_x,...
                indx_interior,indx_edg_W,indx_edg_N,indx_edg_E,indx_edg_S,...
                indx_cor_WN,indx_cor_EN,indx_cor_ES,indx_cor_WS, ...
                indx_mat,ny,nx,ncells,u,alpha_uv);

    end

    
    % #3.- SOLVE  Y MOMENTUM
    A_Ptempv=A_P;% temporal central coeffitients matrix for v

    A_Ptempv(ny,nxSymcond)=A_Pv(nxSymcond);%Replace coeffitients that 
    % should not be the same for the Simetryc boundary condition

    switch solver_y_momentum
       case "G-S"

           [v,rsid_y,err_y,it_inny] = GaussSeidel_GeneralSolver(A_Ptempv,A_W,...
               A_N,A_E,A_S,Su_y,rsid_y,epsilon_v,err_y,max_iterations_v,...
               alpha_uv,nx,ny,v,v);

       case "GMRES"
           [v,rsid_y,err_y,err_y2,it_inny] = gmresSolver(A_Ptempv,A_W,...
               A_N,A_E,A_S,Su_y,A_ysp,suY_vec,restart_y,tol_y,maxit_y,...
               indx_interior,indx_edg_W,indx_edg_N,indx_edg_E,indx_edg_S,...
               indx_cor_WN,indx_cor_EN,indx_cor_ES,indx_cor_WS, ...
               indx_mat,ny,nx,ncells,v,alpha_uv);

       case "BICGSTAB"
           [v,rsid_y,err_y,err_y2,it_inny] = bicgstabSolver(A_Ptempv,A_W,...
               A_N,A_E,A_S,Su_y,A_ysp,suY_vec,restart_y,tol_y,maxit_y,...
               indx_interior,indx_edg_W,indx_edg_N,indx_edg_E,indx_edg_S,...
               indx_cor_WN,indx_cor_EN,indx_cor_ES,indx_cor_WS, ...
               indx_mat,ny,nx,ncells,v,alpha_uv);
    end
    
    %Apply underelaxation to main diagonal coeffitients
    A_P=A_P/alpha_uv;
    A_Pv=A_Pv/alpha_uv;

    % #4.- FACE VELOCITY COMPUTATION USING RHIE-CHOW INTERPOLATION

    [u_f,v_f]=face_vel_intRC(u,v,p,uVecNormFaces,distNeighbNods,...
        uVecNeighbNods,weightDistFactors,grad_p,A_P,A_Pv,cellVols,nx,...
        ny,solidMask);

    %Apply again inlet condition
    u_f(:,1)=u0;

    %Apply Neumman to East and North Edges
    u_f(:,end)=u(:,end);
    v_f(1,:)=-v(1,:);

    % #5.- PRESSURE CORRECTION EQUATION LINK COEFFITIENTS
    [Ap_W,Ap_N,Ap_E,Ap_S,Ap_P,Su_p]=pressureCorr_link_coeff(u_f,v_f,A_P,...
        A_Pv,weightDistFactors,cellVols,distNeighbNods,lgtFaces,nx,ny,...
        solidMask);

    % #6.- SOLVE PRESSURE CORRECTION EQUATION

    %restart pressure correction
    p_prime=zeros(ny,nx);

    switch solver_p_correction

        case "G-S"
            [p_prime,rsid_p,err_p,int_innp]=GaussSeidel_GeneralSolver(...
                Ap_P,Ap_W,Ap_N,Ap_E,Ap_S,Su_p,rsid_p,epsilon_p,err_p,...
                max_iterations_p,1,nx,ny,p_prime,p_prime);

        case "GMRES"
            [p_prime,rsid_p,err_p,err_p2,int_innp] = gmresSolver(Ap_P,...
                Ap_W,Ap_N,Ap_E,Ap_S,Su_p,A_psp,suP_vec,restart_p,tol_p,...
                maxit_p,indx_interior,indx_edg_W,indx_edg_N,indx_edg_E,...
                indx_edg_S,indx_cor_WN,indx_cor_EN,indx_cor_ES,...
                indx_cor_WS,indx_mat,ny,nx,ncells,p_prime,1);

        case "BICGSTAB"
            [p_prime,rsid_p,err_p,err_p2,int_innp] = bicgstabSolver(Ap_P,...
                Ap_W,Ap_N,Ap_E,Ap_S,Su_p,A_psp,suP_vec,restart_p,tol_p,...
                maxit_p,indx_interior,indx_edg_W,indx_edg_N,indx_edg_E,...
                indx_edg_S,indx_cor_WN,indx_cor_EN,indx_cor_ES,...
                indx_cor_WS,indx_mat,ny,nx,ncells,p_prime,1);


    end

    % #7.- CORRECT PRESSURE
    p =p + alpha_p*p_prime;

     % #9.- CORRECT CELL CENTER VELOCITY
    %Gradient of pressure correction field
    [grad_p_prime] = computepPrimeGradient(p_prime,wlsqOperator,nx,ny);

    [u,v]=cvel_correct(A_P,A_Pv,u,v,grad_p_prime,cellVols,1,nx,ny,solidMask);

    % #8.- CORRECT FACE VELOCITY
    [u_f,v_f]=fvel_correct(u_f,v_f,p_prime,A_P,A_Pv,cellVols,...
        weightDistFactors,distNeighbNods,1,nx,ny,solidMask);

    %---------------------  II TURBULENCE   ------------------------------

    %Compute velocity gradients
    [grad_u] = computeuVelGradient(u,wlsqOperator,u0,solidMask,nx,ny);
    [grad_v] = computevVelGradient(v,wlsqOperator,nx,ny);

    %{} this block from here to deactivate turbulence computation
   
    %Compute vorticity*
    [vorticity] = computeVorticity(grad_u,grad_v,ny,nx);
    lower_omega=1.4142*(abs(vorticity));

    %Compute nu_tilde Gradient

    [grad_nu_tilde] = computeNuTGradient(nu_tilde,wlsqOperator,...
          nu,nx,ny,solidMask);

    %Compute nu_tilde at cell vertexes

    nu_tilde_cor=computePhiVertex(nu_tilde,grad_nu_tilde,VecCentVertx,nx,ny);

    %10.- SPALART ALLMARAS MODEL LINK COEFFITIENTS AND SOURCES

    [Ant_W,Ant_N,Ant_E,Ant_S,Ant_P,Su_nt] = saTurbulence_link_coeff( ...
        nu_tilde,nu,lower_omega,lgtFaces,cellVols,u_f,v_f,distMinWall, ...
        grad_nu_tilde,dW,dN,dE,dS,dW_c,dN_c,dE_c,dS_c,nu_tilde_cor,kappa,...
        sigma_sa,cb1,cb2,cv1,cw1,cw2,cw3,nx,ny,solidMask);

    %11.- SOLVE NU_TIDE

    switch solver_nu_tilde

        case "G-S"
            [nu_tilde,rsid_nt,err_nt,int_innnt]=GaussSeidel_GeneralSolver(...
                Ant_P,Ant_W,Ant_N,Ant_E,Ant_S,Su_nt,rsid_nt,epsilon_nt,err_nt,...
                max_iterations_nt,1,nx,ny,nu_tilde,nu_tilde);

        case "GMRES"
            [nu_tilde,rsid_nt,err_nt,err_nt2,int_innnt] = gmresSolver(Ant_P,...
                Ant_W,Ant_N,Ant_E,Ant_S,Su_nt,A_ntsp,suNt_vec,restart_nt,tol_nt,...
                maxit_nt,indx_interior,indx_edg_W,indx_edg_N,indx_edg_E,...
                indx_edg_S,indx_cor_WN,indx_cor_EN,indx_cor_ES,...
                indx_cor_WS,indx_mat,ny,nx,ncells,p_prime,1);

        case "BICGSTAB"
            [nu_tilde,rsid_nt,err_nt,err_nt2,int_innnt] = bicgstabSolver(Ant_P,...
                Ant_W,Ant_N,Ant_E,Ant_S,Su_nt,A_ntsp,suNt_vec,restart_nt,tol_nt,...
                maxit_nt,indx_interior,indx_edg_W,indx_edg_N,indx_edg_E,...
                indx_edg_S,indx_cor_WN,indx_cor_EN,indx_cor_ES,...
                indx_cor_WS,indx_mat,ny,nx,ncells,p_prime,1);
    end


    %12.- COMPUTE EDDY VISCOSITY FROM NU_TILDE

    [mu_turbulent] = saEddyViscosity(nu_tilde,nu,rho,cv1);


    %Store Reynolds,Max-Vorticity, and Max-Eddy viscosity to monitore 
    %turbulence model
    max_abs_vort=max(abs(vorticity(:)));
    max_eddy_visc=max(abs(mu_turbulent(:)));

    turbulenceMonitoring(iterations_cont,:)=[Re,max_abs_vort,...
        max_eddy_visc];

    %}

    %---------------------    III RAMPING   ------------------------------
    
    
    if iterations_cont>=iterations_trig
        %(iterations_cont>=iterations_trig && iterations_cont <20000) || (iterations_cont>30000) % Trigger new condition at certain iteration
        
        
        if ramping_cont<ramping_steps%Stop rampig crteria
            if mod(iterations_cont,ramping_spacing) == 0
                    ramping_cont=ramping_cont+1;
                    % Do something special every n iterations
                    u0=u0*f_ramp;% Ramping U0
                    %enforce again boundary inlet condition
                    u_face(:,1)=u0;
                    Re=u0*rho*(bumpLgt-0.5)/mu; % RECOMPUTE RE
            end
        end

        %}

        
        if mod(iterations_cont,5000) == 0
            file_var_name="ReSolveBoundTurbItNum_" + string(iterations_cont) + ".mat";
            save(file_var_name)
        end
        
    end
    %}
    

    
  
    
    %-----------------  IV EVALUATE CONVERGENCE CRITERIA  -----------------

    %check continuity
    [rsid_cont,err_cont] = checkContinuity(u_f,v_f,lgtFaces,nx,ny);
    
    %Store residuals data
    residualsMat(iterations_cont,:)=[err_x,err_y,err_p,err_nt,err_cont];
    inneritcontMat(iterations_cont,:)=[it_innx,it_inny,int_innp,int_innnt];
    gmrsErrors(iterations_cont,:)=[err_x2,err_y2,err_p2,err_nt2];
    disp(residualsMat(iterations_cont,:)) %display current error

    if err_x < error_tgt && err_y < error_tgt && err_cont < error_tgt ...
            && err_p < error_tgt
       convergedFlg=true;
       fprintf("Converged  at iteration\n")
       disp(iterations_cont)
       
    elseif err_x>max_residual || err_y > max_residual ||...
            err_cont > max_residual || err_p>max_residual
        fprintf("Unstable solution iterations stopped at iteration ")
        disp(iterations_cont)
        break
    elseif iterations_cont >= max_iterations
        fprintf("Max iterations reached \n")
        disp(iterations_cont)
        break
    elseif isnan(err_x) || isnan(err_y) || isnan(err_cont) || isnan(err_p)
        fprintf("The system became undetermined  at iteration \n")
        disp(iterations_cont)
        break
    else
        convergedFlg=false;
    end

    %convergedFlg=true;

end
toc
%end loop

%% POST PROCESS

%VELOCITY MAGNITUDE
v_m= sqrt(u.^2 + v.^2);

fprintf("Max Vorticity  =  ")
disp(max(abs(vorticity(:))));

%PRESSURE COEFFITIENT
x_pressure_coeff=zeros(size(nxSolid));

for j=1:nx_fp
    j_1=nxSolid(j);
    x_pressure_coeff(j)=faceCentrs(ny,j_1,4,1);
end 
x_pressure_coeff=x_pressure_coeff-(domLgt-bumpLgt)/2;
pressure_coeff=(p(ny,nxSolid)-p0)/(0.5*rho*u0^2);

%Compute gradient of Reynolds stresses
grad_tau_xx=computeTauGradient(tau_xx,wlsqOperator,grad_tau_xx);
grad_tau_xy=computeTauGradient(tau_xy,wlsqOperator,grad_tau_xy);
grad_tau_yy=computeTauGradient(tau_yy,wlsqOperator,grad_tau_yy);

%compute turbulent stresses
[tau_xx,tau_xy,tau_yy] = saBoussinesqTurbulentStresses(...
    mu_turbulent,grad_u,grad_v,tau_xx,tau_xy,tau_yy,nx,ny);

%Nondimensional Eddy viscosity
nondim_mu_turbulent=mu_turbulent/mu;

%% PLOTS
%RESIDUALS 

figure(2)
plot(1:iterations_cont-1,...
    residualsMat(1:iterations_cont-1,1),...
    1:iterations_cont-1,...
    residualsMat(1:iterations_cont-1,2),...
    1:iterations_cont-1,...
    residualsMat(1:iterations_cont-1,3),...
    1:iterations_cont-1,...
    residualsMat(1:iterations_cont-1,4),...
    1:iterations_cont-1,...
    residualsMat(1:iterations_cont-1,5))
legend("u vel","v vel","Pressure correction","Turbulence","Continuity")
title("Residuals")
xlabel("Iterations")
ylabel("Residual")
yscale log

figure(3) %INNER ITERATIONS
plot(1:iterations_cont-1,...
    inneritcontMat(1:iterations_cont-1,1),...
    1:iterations_cont-1,...
    inneritcontMat(1:iterations_cont-1,2),...
    1:iterations_cont-1,...
    inneritcontMat(1:iterations_cont-1,3),...
    1:iterations_cont-1,...
    inneritcontMat(1:iterations_cont-1,4))
legend("u vel","v vel","Pressure correction","turbulence")
title("Inner iterations per outer iteration")
xlabel("Iterations")
ylabel("Inner iterations")

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

figure(6) %VELOCITY MAGNITUDE
contourf(Xctrs,Yctrs,v_m, 20, 'LineColor', 'none')
title("Velocity Magnitude (m/s)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

figure(7) %pressure correction
contourf(Xctrs,Yctrs,p_prime, 20, 'LineColor', 'none')
title("Pressure correction (N/m2)")
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

figure(9) %Vorticity
contourf(Xctrs,Yctrs,vorticity, 20, 'LineColor', 'none')
title("Vorticity (1/s)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

figure(10)%Streamlines
streamslice(Xctrs, Yctrs, u, v); % For a sparse representation
title("streamlines")
xlabel("Lenght (m)")
ylabel("Height (m)")
axis equal


figure(11)%Pressure coeffitient
plot(x_pressure_coeff,pressure_coeff ,'- o')
title("PRESSURE COEFFITIENT ON THE PLATE ")
xlabel("X POSITION")
ylabel("PRESSURE COEFFITIENT ")
grid on 
ax = gca;
ax.YDir = 'reverse';


figure(12) %Rotation tensor norm
contourf(Xctrs,Yctrs,lower_omega, 20, 'LineColor', 'none')
title("Rotation tensor norm (1/s)")
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

figure(14)%Eddy viscosity
contourf(Xctrs,Yctrs,mu_turbulent, 20, 'LineColor', 'none')
title("Eddy Viscosity")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

figure(15)%Non Dimensional Eddy viscosity
contourf(Xctrs,Yctrs,mu_turbulent/mu, 20, 'LineColor', 'none')
title("Non Dimensional Eddy Viscosity")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

figure(16) %ITurbulence monitoring RE
plot(1:iterations_cont-1,...
    turbulenceMonitoring(1:iterations_cont-1,1))
title("Reynolds Number")
xlabel("Iterations")
ylabel("Re")

figure(17) %ITurbulence monitoring Max Vorticity
plot(1:iterations_cont-1,...
    turbulenceMonitoring(1:iterations_cont-1,2))
title("Max Vorticity")
xlabel("Iterations")
ylabel("Vorticity (1/s)")

figure(18) %ITurbulence monitoring  Eddy Viscosity
plot(1:iterations_cont-1,...
    turbulenceMonitoring(1:iterations_cont-1,3))
title("Max Eddy viscosity")
xlabel("Iterations")
ylabel("Values")

figure(19)
contourf(Xctrs,Yctrs,tau_xx, 20, 'LineColor', 'none')
title("Turbulent Stress Tau_xx")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

figure(20)
contourf(Xctrs,Yctrs,tau_xy, 20, 'LineColor', 'none')
title("Turbulent Stress Tau_xy")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

figure(21)
contourf(Xctrs,Yctrs,tau_yy, 20, 'LineColor', 'none')
title("Turbulent Stress Tau_yy")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal


%}