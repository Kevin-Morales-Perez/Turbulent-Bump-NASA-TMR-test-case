%RESIDUAL PLOTS
close all

%x-momentum
figure(1)
contourf(Xctrs,Yctrs,rsid_x, 20, 'LineColor', 'none')
title("Raw Residuals for X-momentum")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

%y-momentum
figure(2)
contourf(Xctrs,Yctrs,rsid_y, 20, 'LineColor', 'none')
title("Raw Residuals for Y-momentum")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

%Pressure correction
figure(3)
contourf(Xctrs,Yctrs,rsid_p, 20, 'LineColor', 'none')
title("Raw Residuals for pressure corretion")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

%Continuity
figure(4)
contourf(Xctrs,Yctrs,rsid_cont, 20, 'LineColor', 'none')
title("Continuity Residual")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

%Turbulence
figure(5)
contourf(Xctrs,Yctrs,rsid_nt, 20, 'LineColor', 'none')
title("Turbulence Residual")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal