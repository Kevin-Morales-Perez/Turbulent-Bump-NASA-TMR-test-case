function [u_f,v_f] = face_vel_int(u,v,uVecNormFaces,...
    weightDistFactors,nx,ny)
    %u: x  component velocity field , cell centered 
    %v: y component velocity field, cell centered
    
    %u_face: velocity normal to faces weast and east in the face centers
    %v_face:velocity normal to faces north and south in the face centers
    
    %uVecNormFaces:Array with vectors normal outward to each face of cells
   

    %Interpolation for face velocities 
    %   uf= dot(uP*wfp + uA*(1-wfp),n_f)
    
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
        
            %Apply formula *
            u_f(i,j+1)= u_f_e ;
       
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
            %negative sign to make positive the velocity *
            n_s=reshape(uVecNormFaces(i,j,4,:),[1,2]);

            %get weight factor at cell center for node S *
            weight_S=weightDistFactors(i,j,4);

            %Average velocity vector at face S *
            U_s =weight_S*U_P + (1-weight_S)*U_S;

            %Normal velocity at face S *
            v_f_s=dot(U_s,n_s);
        
            %Apply formula *
            v_f(i+1,j) = v_f_s;

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

        %Apply formula
        v_f(i+1,j) = v_f_s;

    end
end