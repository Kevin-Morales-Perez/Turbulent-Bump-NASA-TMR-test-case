function [s_dc_total_x,s_dc_total_y] = tvdDefCorr(f_C,u_vec,v_vec,...
    grad_u_mat,grad_v_mat,r_mat)
    %TVD deferred corrections
    %Van Leer Limiter
    %Upwind MUSCL
    
    %f_C: Vector with mass fluxes of the cell
    %u_vec,v_vec: Nodal Values 
    %grad_u_mat,grad_v_mat: Arrays with gradients
    %r_mat: Arrays with vectors from cell p to cell A (Neigborhood)
    
   
    %Initialize deferred corrections
    %X
    s_dc_total_x=0;
    %Y
    s_dc_total_y=0;
    
    %Iterate over the array with mass fluxes
    for i=1:length(f_C)
    
        fA=f_C(i);%Mass flux in face A
        
        if fA ~=0
    
            r_P_A=r_mat(i,:);%Get vector from node P  to node A
        
            if fA > 0
        
                r_up_do=r_P_A;%Vector from upwind to downwind node
        
                %X contribution
                phi_x_up=u_vec(end);%upwind node value
                phi_x_do=u_vec(i);%Downwind node value
        
                denom_x=phi_x_do - phi_x_up;%Denominator
    
                if denom_x ~=0%Avoiding division by zero
    
                    %Get upwind cell node gradient
                    grad_x_up=grad_u_mat(end,:);
    
                    %X contribution
                    %Compute Ratio with Darwish formula
                    ratio_x = (2*dot(grad_x_up,r_up_do)/denom_x)-1;

                    %Compute deferred contribution with van leer limiter
                    s_dc_x_a=-0.5*fA*vanLeerLimiter(ratio_x)*denom_x;
                
                    %Sum face contribution to total contribution
                    s_dc_total_x=s_dc_total_x + s_dc_x_a;
    
                end
        
                %Y contribution
                phi_y_up=v_vec(end);%upwind node value
                phi_y_do=v_vec(i);%Downwind node value
    
                denom_y=phi_y_do - phi_y_up;%Denominator
    
                if denom_y ~=0%Avoiding division by zero
    
                    %Get upwind cell node gradient
                    grad_y_up=grad_v_mat(end,:);

                    %Y contribution

                    %Compute Ratio with Darwish formula
                    ratio_y = (2*dot(grad_y_up,r_up_do)/denom_y)-1;

                    %Compute deferred contribution with van leer limiter
                    s_dc_y_a=-0.5*fA*vanLeerLimiter(ratio_y)*denom_y;
                
                    %Sum face contribution to total contribution
                    s_dc_total_y=s_dc_total_y + s_dc_y_a;
    
                end
        
            elseif fA < 0
        
                r_up_do=-r_P_A;%Vector from upwind to downwind node
        
                %X contribution
                phi_x_up=u_vec(i);%upwind node value
                phi_x_do=u_vec(end);%Downwind node value
        
                denom_x=phi_x_do - phi_x_up;%Denominator
    
                if denom_x ~= 0%Avoiding division by zero
    
                    %Get upwind cell node gradient
                    grad_x_up=grad_u_mat(i,:);
    
                    %X contribution
                    %Compute Ratio with Darwish formula
                    ratio_x = (2*dot(grad_x_up,r_up_do)/denom_x)-1;

                    %Compute deferred contribution with van leer limiter
                    s_dc_x_a=-0.5*fA*vanLeerLimiter(ratio_x)*denom_x;

                    %Sum face contribution to total contribution
                    s_dc_total_x=s_dc_total_x + s_dc_x_a;
    
                end
        
                %Y contribution
                phi_y_up=v_vec(i);%upwind node value
                phi_y_do=v_vec(end);%Downwind node value
        
                denom_y=phi_y_do - phi_y_up;%Denominator
    
                if denom_y ~= 0
    
                    %Get upwind cell node gradient
                    grad_y_up=grad_v_mat(i,:);
    
                    %Y contribution

                    %Compute Ratio with Darwish formula
                    ratio_y = (2*dot(grad_y_up,r_up_do)/denom_y)-1;

                    %Compute deferred contribution with van leer limiter
                    s_dc_y_a=-0.5*fA*vanLeerLimiter(ratio_y)*denom_y;
                
                    %Sum face contribution to total contribution
                    s_dc_total_y=s_dc_total_y + s_dc_y_a;
    
                end
        
            end
               
        end   
    
    end

end

