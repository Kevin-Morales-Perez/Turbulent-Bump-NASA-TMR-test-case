%Face Velocity Plots
close all



%Mesh points for u
Xctrs_face_u=zeros(ny,nx+1); 
Yctrs_face_u=zeros(ny,nx+1); 

%Mesh points for v
Xctrs_face_v=zeros(ny+1,nx); 
Yctrs_face_v=zeros(ny+1,nx);

v_face_normal=zeros(ny+1,nx); 
v_face_normal_test=zeros(ny+1,nx);


for i=1:ny
    for j=1:nx

        %U
        Xctrs_face_u(i,j)=faceCentrs(i,j,1,1);
        Xctrs_face_u(i,j+1)=faceCentrs(i,j,3,1);

        Yctrs_face_u(i,j)=faceCentrs(i,j,1,2);
        Yctrs_face_u(i,j+1)=faceCentrs(i,j,3,2);

        %V
        Xctrs_face_v(i,j)=faceCentrs(i,j,2,1);
        Xctrs_face_v(i+1,j)=faceCentrs(i,j,4,1);

        Yctrs_face_v(i,j)=faceCentrs(i,j,2,2);
        Yctrs_face_v(i+1,j)=faceCentrs(i,j,4,2);

        v_face_normal(i,j)=v_f(i,j)*uVecNormFaces(i,j,2,2);
        v_face_normal(i+1,j)=-v_f(i,j)*uVecNormFaces(i,j,4,2);

        v_face_normal_test(i,j)=v_face_test(i,j)*uVecNormFaces(i,j,2,2);
        v_face_normal_test(i+1,j)=-v_face_test(i,j)*uVecNormFaces(i,j,4,2);


    end
end



figure(1) %CELL FACE VELOCITY IN X AXIS
contourf(Xctrs_face_u,Yctrs_face_u,u_f, 20, 'LineColor', 'none')
title("Face Velocity in x axis (m/s)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

figure(2) %CELL  FACE VELOCITY IN Y AXIS
contourf(Xctrs_face_v,Yctrs_face_v,v_face_normal, 20, 'LineColor', 'none')
title("Velocity in y axis (m/s)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

figure(3) %CELL  FACE VELOCITY IN Y AXIS
contourf(Xctrs_face_v,Yctrs_face_v,v_face_normal_test, 20, 'LineColor', 'none')
title("Velocity in y axis Test (m/s)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal