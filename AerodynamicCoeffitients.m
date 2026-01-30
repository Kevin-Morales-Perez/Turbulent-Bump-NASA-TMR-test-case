% Aerodynamic Coeffitients 



%% Surface skin friction coeffitient
x_surfSkinFricCoeff=zeros(size(nxSolid));%X coordinates 

for j=1:nx_fp
    j_1=nxSolid(j);
    x_surfSkinFricCoeff(j)=faceCentrs(ny,j_1,4,1);
end 

x_surfSkinFricCoeff=x_surfSkinFricCoeff-(domLgt-bumpLgt)/2;%Adjust to grid
%global coordinates

%Velocity in cells adyacent to surface
u_surface=u(ny,nxSolid);
v_surface=v(ny,nxSolid);

%distance from cell center to surface
d_surface=distNeighbNods(ny,nxSolid,4);%distMinWall(ny,nxSolid);

%Effective viscosity
mu_effective_surface=mu_turbulent(ny,nxSolid) + mu;

%Velocity tangent to south cell face 
u_surface_norm=zeros(size(nxSolid));

%dynamic pressure
dynPress=0.5*rho*u0^2;

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


surfSkinFricCoeff=shear_surface_stress/dynPress;

%pressure_coeff=(p(ny,nxSolid)-p0)/(0.5*rho*u0^2);

figure(1)%Surface Skin Friction Coeffitient 
plot(x_surfSkinFricCoeff , surfSkinFricCoeff ,'- o')
title("SURFACE SKIN FRICTION COEFFITIENT")
xlabel("X POSITION")
ylabel("Cf")
grid on 
axis([0 1.5 0 0.008])
%}

%% PRESSURE COEFFITIENT
x_pressure_coeff=zeros(size(nxSolid));%X coordinates for pressure coeffitient

for j=1:nx_fp
    j_1=nxSolid(j);
    x_pressure_coeff(j)=faceCentrs(ny,j_1,4,1);
end 
x_pressure_coeff=x_pressure_coeff-(domLgt-bumpLgt)/2;


pressure_coeff=(p(ny,nxSolid)-p0)/(0.5*rho*u0^2);