function [u_f,v_f] = face_vel_intRC(u,v,p,uVecNormFaces,distNeighbNods,...
    uVecNeighbNods,weightDistFactors,grad_p,aP,aPv,cellVols,nx,ny,solidMask)
    %u: x  component velocity field , cell centered 
    %v: y component velocity field, cell centered
    %p: pressure cell centered
    %u_face: velocity normal to faces weast and east in the face centers
    %v_face:velocity normal to faces north and south in the face centers
    %VecCentNodFaceCents: Array with vectors from cell center to face
    %centers
    %uVecNormFaces:Array with vectors normal outward to each face of cells
    %distNeighbNods: Distances between cell centers and neigborhood nodes
    %uVecNeighbNods:Unitary vectors from cell center to nb cell centers
    %grad_u: Gradient for u 
    %grad_v: Gradient for v
    %aP:Coeffitient from momentum equations
    %aPv:Coeffitient from momentum equations for symetric Boundary
    %condition
    %cellVols: volumes of the cells
    %nx_upstr: number of cells in the free upwind zone
    %nx_dwnstr: number of cells in the free downstream zone 

        %Rie chow interpolation for face velocities 
        %   uf= dot(uP*wfp + uA*(1-wfp),n_f) + df*(Pp-Pa) - 0.5dot(((VP/aP)*grad(p)p + 
        %  (VA/aA)*grad(p)A)),epsilon_f)
    %Orthogonal - and non orthogonal meshes

    %Initalize variables
    u_f=zeros(ny,nx+1);
    v_f=zeros(ny+1,nx);

    % west and east face velocities  (u_face) *
    
    for j=1:nx-1
        for i=1:ny
        
            %Compute average velocity_____________________________________
                      
            %get both velocity components at the cell center P *       
            u_p=u(i,j);
            v_p=v(i,j);
            %group in a vector
            U_P=[u_p,v_p];
            %get both velocity components at the next cell center E*
            u_e=u(i,j+1);       
            v_e=v(i,j+1);

            %group in a vector *
            U_E=[u_e,v_e];
      
            %get unitary vector normal to face east *        
            n_e=reshape(uVecNormFaces(i,j,3,:),[1,2]);
        
            %get weight factor at cell center for node E *
            weight_E=weightDistFactors(i,j,3);

            %Averaged velocity vector at face E * 
            U_e=weight_E*U_P + (1-weight_E)*U_E;

            %Normal velocity at face E* 
            u_f_e=dot(U_e,n_e);
        
            %use Rie Chow interpolation formula___________________________
            
            %get distance from cell center to east cell center *
            dist_p_e =distNeighbNods(i,j,3);
        
            %get unitary vector from cell center to east nb cell center  *       
            uvec_p_e=reshape(uVecNeighbNods(i,j,3,:),[1,2]);
          
            %Pressure gradients*
            gradpP_x=grad_p(i,j,1,1);
            gradpP_y=grad_p(i,j,1,2);
            gradpP=[gradpP_x,gradpP_y];
            %________________________*
            gradpE_x=grad_p(i,j+1,1,1);
            gradpE_y=grad_p(i,j+1,1,2);
            gradpE=[gradpE_x,gradpE_y];

            %Averaged Coeffitient for pressure gradient at face e*
            De_av=(weight_E*(cellVols(i,j)/aP(i,j)) +...
                (1 -weight_E)*(cellVols(i,j+1)/aP(i,j+1)));

            %Average pressure gradient at face e *
            gradpe_av=(cellVols(i,j)/aP(i,j))*gradpP*weight_E + ...
                (cellVols(i,j+1)/aP(i,j+1))*gradpE*(1 -weight_E);

            %Apply formula *
            u_f(i,j+1)= u_f_e + De_av*((p(i,j)-p(i,j+1))/dist_p_e) - ...
                dot(gradpe_av,uvec_p_e);
       
        end
    end

    % North and South face velocities  (v_face)
    
    for i=1:ny-2
        for j=1:nx

            %Compute average velocity____________________________________
        
            %get both velocity components at the cell center P *       
            u_p=u(i,j);
            v_p=v(i,j);
            %group in a vector
            U_P=[u_p,v_p];
            %get both velocity components at the next cell center S *
            u_s=u(i+1,j);        
            v_s=v(i+1,j);

            %group in a vector *
            U_S=[u_s,v_s];

            %get unitary vector normal to face south *
            n_s=reshape(uVecNormFaces(i,j,4,:),[1,2]);

            %get weight factor at cell center for node S *
            weight_S=weightDistFactors(i,j,4);

            %Average velocity vector at face S *
            U_s =weight_S*U_P + (1-weight_S)*U_S;

            %Normal velocity at face S *
            v_f_s=dot(U_s,n_s);

            %use Rie Chow interpolation formula___________________________

            %get distance from cell center to south cell center *
            dist_p_s = distNeighbNods(i,j,4);
        
            %get normal vector from cell center to south cell center *
            uvec_p_s = reshape(uVecNeighbNods(i,j,4,:), [1, 2]);

            %Pressure gradients *
            gradpP_x=grad_p(i,j,1,1);
            gradpP_y=grad_p(i,j,1,2);
            gradpP=[gradpP_x,gradpP_y];
            %__________________________
            gradpS_x=grad_p(i+1,j,1,1);
            gradpS_y=grad_p(i+1,j,1,2);
            gradpS = [gradpS_x,gradpS_y];

            %Coeffitient for pressure gradient at face s
            Ds_av=(weight_S*(cellVols(i,j)/aP(i,j)) +...
                (1 -weight_S)*(cellVols(i+1,j)/aP(i+1,j)));

            %Average pressure gradient at face s *
            gradps_av=(cellVols(i,j)/aP(i,j))*gradpP*weight_S + ...
                (cellVols(i+1,j)/aP(i+1,j))*gradpS*(1 -weight_S);
        
            %Apply formula *
            v_f(i+1,j) = v_f_s + Ds_av*((p(i,j)-p(i+1,j))/dist_p_s) - ...
                dot(gradps_av,uvec_p_s);

        end
    end

    %South Edge*
    i=ny-1;

    for j=1:nx
        %Compute average velocity____________________________________
        
        %get both velocity components at the cell center P        
        u_p=u(i,j);
        v_p=v(i,j);
        %group in a vector
        U_P=[u_p,v_p];
        %get both velocity components at the next cell center S
        u_s=u(i+1,j);        
        v_s=v(i+1,j);

        %group in a vector
        U_S=[u_s,v_s];

        %get unitary vector normal to face south
        %negative sign to make positive the velocity
        n_s=reshape(uVecNormFaces(i,j,4,:),[1,2]);

        %get weight factor at cell center for node S
        weight_S=weightDistFactors(i,j,4);

        %Average velocity vector at face S
        U_s =weight_S*U_P + (1-weight_S)*U_S;

        %Normal velocity at face S
        v_f_s=dot(U_s,n_s);

        %use Rie Chow interpolation formula___________________________

        %get distance from cell center to south cell center
        dist_p_s = distNeighbNods(i,j,4);
    
        %get normal vector from cell center to south cell center 
        uvec_p_s = reshape(uVecNeighbNods(i,j,4,:), [1, 2]);

        %Pressure gradients
        gradpP_x=grad_p(i,j,1,1);
        gradpP_y=grad_p(i,j,1,2);
        gradpP=[gradpP_x,gradpP_y];
        %__________________________
        gradpS_x=grad_p(i+1,j,1,1);
        gradpS_y=grad_p(i+1,j,1,2);
        gradpS = [gradpS_x,gradpS_y];

        if solidMask(j) 
            
            %Normal behaviour (Solid surface)

            %Coeffitient for pressure gradient at face s
            Ds_av=weight_S*(cellVols(i,j)/aP(i,j)) +...
                (1 -weight_S)*(cellVols(i+1,j)/aP(i+1,j));

            %Average pressure gradient at face s
            gradps_av=(cellVols(i,j)/aP(i,j))*gradpP*weight_S + ...
                (cellVols(i+1,j)/aP(i+1,j))*gradpS*(1 -weight_S);
    
            
        else

            %use aPv coeffitients (Free stream zones)
            %Coeffitient for pressure gradient at face s
            Ds_av=weight_S*(cellVols(i,j)/aP(i,j)) +...
                (1 -weight_S)*(cellVols(i+1,j)/aPv(j));

            %Average pressure gradient at face s
            gradps_av=(cellVols(i,j)/aP(i,j))*gradpP*weight_S + ...
                (cellVols(i+1,j)/aPv(j))*gradpS*(1 -weight_S);
            
        end

        %Apply formula
        v_f(i+1,j) = v_f_s + Ds_av*((p(i,j)-p(i+1,j))/dist_p_s) - ...
            dot(gradps_av,uvec_p_s);

    end
end
