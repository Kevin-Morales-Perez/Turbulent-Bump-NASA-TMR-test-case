function [vorticity] = computeVorticity(grad_u,grad_v,ny,nx)
    
%Initialize variable
vorticity=zeros(ny,nx);


    %Computation of vorticity
    for i=1:ny
        for j=1:nx
            dv_dx=grad_v(i,j,1,1);
            du_dy=grad_u(i,j,1,2);
            vorticity(i,j)=dv_dx-du_dy;           
        end 
    end
    
end