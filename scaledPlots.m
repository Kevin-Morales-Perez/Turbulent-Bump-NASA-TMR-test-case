%TURBULENT PLOTS SCALED

close all

figure(1) %VELOCITY MAGNITUDE
contourf(Xctrs,Yctrs,v_m, 20, 'LineColor', 'none')
title("Velocity Magnitude (m/s)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis([distPlat distPlat+1.5 0 0.08])

figure(2)%Nu tilde
contourf(Xctrs,Yctrs,nu_tilde, 20, 'LineColor', 'none')
title("Spalart - Allmaras nu~ transported variable")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis([distPlat distPlat+1.5 0 0.08])

figure(3)%Eddy viscosity
contourf(Xctrs,Yctrs,mu_turbulent, 20, 'LineColor', 'none')
title("Eddy Viscosity")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis([distPlat distPlat+1.5 0 0.08])


figure(4)%Rotation tensor norm
contourf(Xctrs,Yctrs,lower_omega, 20, 'LineColor', 'none')
title("Rotation tensor norm")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis([distPlat distPlat+1.5 0 0.08])
   
figure(5)%Non dimensional Eddy viscosity
contourf(Xctrs,Yctrs,nondim_mu_turbulent, 20, 'LineColor', 'none')
title("Non dimensional eddy viscosity")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis([distPlat distPlat+1.5 0 0.08])
