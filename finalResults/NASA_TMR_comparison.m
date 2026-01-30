%Comparison between results of this code and NASA turbulence modelling
%resource 

%This code  VS FUN3D (Fully Unstructured Navier - Stokes 3D)
%8500 vs 903169 cells

%( FUN3D is a Computational Fluid Dynamics (CFD) suite of tools actively
%  developed at NASA)
close all
%Load converged result from this code
%load("ReSolveBoundTurbItNum_900000.mat")
load("finalConvergedRe3e6.mat")
%load("Fixing2TurbItNum_280000.mat")
%% PRESSURE COEFITIENT 

q=0.5*rho*u0^2;%dynamic pressure 

%This code 
x_surface_coeff=zeros(size(nxSolid));%X coordinates for pressure coeffitient

for j=1:nx_fp
    j_1=nxSolid(j);
    x_surface_coeff(j)=faceCentrs(ny,j_1,4,1);
end 
x_surface_coeff=x_surface_coeff-(domLgt-bumpLgt)/2;%Adjust to grid 


pressure_coeff=(p(ny,nxSolid)-p0)/(q);

%FUN3D

load("fun3dcp.mat")

figure(1)%Pressure coeffitient
plot(x_surface_coeff,pressure_coeff ,'- o',fun3d_cp.X,fun3d_cp.CP ,'- *')
legend({'Código propio', 'NASA CFL3D/FUN3D'}, 'Location', 'best')
title('Coeficiente de presión sobre la protuberancia (Re = 3\times10^6)')
xlabel('Posición en x (m)')
ylabel('C_p')
grid on 
ax = gca;
ax.YDir = 'reverse';


%% Surface skin friction coeffitient

%Velocity in cells adyacent to surface
u_surface=u(ny,nxSolid);
v_surface=v(ny,nxSolid);

%distance from cell center to surface
d_surface=distNeighbNods(ny,nxSolid,4);%distMinWall(ny,nxSolid);

%Effective viscosity
mu_effective_surface=mu_turbulent(ny,nxSolid) + mu;

%Velocity tangent to south cell face 
u_surface_norm=zeros(size(nxSolid));

%Computing tangential velocity to surface 
for j=1:nx_fp
    j_2=nxSolid(j);

    %Unitary Vector tangential to the surface 
    v_t_surf=reshape(uVecParlFaces(ny,j_2,4,:),[1,2]);

    %velocity vector
    u_vctr =[u_surface(j),v_surface(j)];

    u_surface_norm(j)=dot(v_t_surf,u_vctr);
   
end

shear_surface_stress=mu_effective_surface.*(u_surface_norm./d_surface);

surfSkinFricCoeff=shear_surface_stress/q;


%load
load("fun3dcf.mat")

figure(2)%Surface Skin Friction Coeffitient 
plot(x_surface_coeff , surfSkinFricCoeff ,'- o',fun3d_cf.X,fun3d_cf.CF ,'- *')
legend('Código propio', 'NASA CFL3D/FUN3D', 'Location', 'best')
title('Coeficiente de fricción sobre la protuberancia')
xlabel('Posición en x')
ylabel('C_f')
grid on 
axis([0 1.5 0 0.008])

%% LIFT AND DRAG
%nx_fp: points in the solid plate
%nxSolid: indexes for cells above solid plate

%Direction of velocity
u_direction = [1,0];

y_positive_direction=[0,1];

pressure_lift_vec =zeros(1,nx_fp);
pressure_drag_vec= zeros(1,nx_fp);

friction_lift_vec =zeros(1,nx_fp);
friction_drag_vec= zeros(1,nx_fp);

shear_lift_vec=zeros(1,nx_fp);



for j = 1:nx_fp 
    j_1=nxSolid(j);

    %Pressure forces
    %get pressure coeffitient
    cpi=pressure_coeff(j);
    
    %get vector normal to face S (Face south of cell above solid plate)
    nsi = reshape(uVecNormFaces(ny,j_1,4,:),[1,2]);
    
    %get face area
    fsi=lgtFaces(ny,j_1,4);

    %Normal Force in the cell face
    wsi =cpi*fsi*nsi;

    %lift and drag
    pressure_lift_vec(j)=dot(wsi,y_positive_direction);
    pressure_drag_vec(j)=dot(wsi,u_direction);


    %Shear forces
    %get skin friction coeffitient
    cfi=surfSkinFricCoeff(j);

    %get vector paralel to face S
    esi = reshape(uVecParlFaces(ny,j_1,4,:),[1,2]);

    %shear force in the cell face
    tsi=cfi*fsi*esi;

    %lift and drag
    friction_lift_vec(j)=dot(tsi,y_positive_direction);
    friction_drag_vec(j)=dot(tsi,u_direction);

end

pressure_lift =q*sum(pressure_lift_vec);
pressure_drag =q*sum(pressure_drag_vec);

cl_p=pressure_lift/(q*bumpLgt);%Pressure Lift coefficient 
cd_p=pressure_drag/(q*bumpLgt);%Pressure Drag Coefficient

shear_lift =q*sum(friction_lift_vec);
shear_drag =q*sum(friction_drag_vec);

cl_v=shear_lift/(q*bumpLgt);%Friction Lift coefficient 
cd_v=shear_drag/(q*bumpLgt);%Friction Drag Coefficient


cl=cl_p+cl_v;
cd=cd_p+cd_v;


%% PRESSURE TOTAL FORCE

Fp=zeros(nx_fp,2);

for j = 1:nx_fp 
    j_1=nxSolid(j);

    %get pressure coeffitient
    cpi=pressure_coeff(j);

    %get face area
    fsi=lgtFaces(ny,j_1,4);

    %get vector normal to face S (Face south of cell above solid plate)
    nsi = reshape(uVecNormFaces(ny,j_1,4,:),[1,2]);

    Fp(j,:)=cpi*q*fsi*nsi;

end

plot(x_surface_coeff,Fp(:,2));

Fp_l=sum(Fp(:,2));
Fp_d=sum(Fp(:,1));