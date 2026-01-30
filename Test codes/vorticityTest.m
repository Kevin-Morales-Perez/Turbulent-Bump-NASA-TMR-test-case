%% ADITIONAL TEMPORAL PLOTS

grad_u_x=zeros(ny,nx);
grad_u_y=zeros(ny,nx);
grad_v_x=zeros(ny,nx);
grad_v_y=zeros(ny,nx);

for i=1:ny
    for j=1:nx

        grad_u_x(i,j)=grad_u(i,j,1,1);
        grad_u_y(i,j)=grad_u(i,j,1,2);
        grad_v_x(i,j)=grad_v(i,j,1,1);
        grad_v_y(i,j)=grad_v(i,j,1,2);


    end
end



figure(11) %Grad_u_x
contourf(Xctrs,Yctrs,grad_u_x, 20, 'LineColor', 'none')
title("grad_u_x (N/m2)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal


figure(12) %Grad_u_x
contourf(Xctrs,Yctrs,grad_u_y, 20, 'LineColor', 'none')
title("grad_u_y (N/m2)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal


figure(13) %Grad_u_x
contourf(Xctrs,Yctrs,grad_v_x, 20, 'LineColor', 'none')
title("grad_v_x (N/m2)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

figure(14) %Grad_u_x
contourf(Xctrs,Yctrs,grad_v_y, 20, 'LineColor', 'none')
title("grad_v_y (N/m2)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

figure(15) %vorticity
contourf(Xctrs,Yctrs,vorticity, 20, 'LineColor', 'none')
title("grad_v_y (N/m2)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal