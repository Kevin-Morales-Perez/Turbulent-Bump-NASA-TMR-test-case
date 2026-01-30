%Code to pre process interpolated fields to main sover

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


% 0.1- compute pressure gradient *
grad_p=computePressGradient(p,wlsqOperator,p0,f1p,f2p,nx,ny);
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
%12.- COMPUTE EDDY VISCOSITY FROM NU_TILDE

[mu_turbulent] = saEddyViscosity(nu_tilde,nu,rho,cv1);
% #4.- FACE VELOCITY COMPUTATION 

[u_f,v_f] = face_vel_int(u,v,uVecNormFaces,...
    weightDistFactors,nx,ny);

%Apply again inlet condition
u_f(:,1)=u0;

%Apply Neumman to East and North Edges
u_f(:,end)=u(:,end);
v_f(1,:)=-v(1,:);
%check continuity
[rsid_cont,err_cont] = checkContinuity(u_f,v_f,lgtFaces,nx,ny);


save("NewInterpolationtomesh13.mat")
