function [u,v] = cvel_correct(aP,aPv,u_star,v_star,grad_p_prime,...
    cellVols,alpha_uv,nx,ny,solidMask)
%CORRECTION OF CELL CENTER VELOCITIES 
% using pressure corection field to correct CELL CENTER velocities
%u_star, v_star: velocities from previous iteration
%aP, aPv:Coeffitient from momentum equations  
%u,v: velocities at nodes
%grad_p_prime: gradient for pressure correction
%cellVols: volumes of the cells

    u=u_star;
    v=v_star;
    
    %Velocity in X and Y axis *
    for i = 1:ny-1
        for j = 1:nx
            
            %Get pressure correction gradient components*
            gradP_px=grad_p_prime(i,j,1,1);
            gradP_py=grad_p_prime(i,j,1,2);
            
            u(i,j)=u_star(i,j) - alpha_uv*(cellVols(i,j)/aP(i,j))*gradP_px;
            v(i,j)=v_star(i,j) - alpha_uv*(cellVols(i,j)/aP(i,j))*gradP_py;
       
        end
    end

    i=ny;
    
    for j=1:nx

        %Get pressure correction gradient components*
        gradP_px=grad_p_prime(i,j,1,1);
        gradP_py=grad_p_prime(i,j,1,2);

        u(i,j)=u_star(i,j) - alpha_uv*(cellVols(i,j)/aP(i,j))*gradP_px;
        
        if solidMask(j)%Solid Surface
       
            v(i,j)=v_star(i,j) - alpha_uv*(cellVols(i,j)/aP(i,j))*gradP_py;
            
        else%Symetric condition

            v(i,j)=v_star(i,j) - alpha_uv*(cellVols(i,j)/aPv(j))*gradP_py;
            
        end
    end

