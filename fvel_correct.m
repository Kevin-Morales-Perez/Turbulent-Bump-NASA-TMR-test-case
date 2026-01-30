function [u_f,v_f] = fvel_correct(u_f_star,v_f_star,p_prime,aP,aPv,...
    cellVols,weightDistFactors,distNeighbNods,alpha_uv,nx,ny,solidMask)
%FACE VELOCITY CORRECTION
%u_face,v_face: velocities normal to faces
%p_prime: pressure correction field at cell centers
%%aP, aPv:central Coeffitient from momentum equations  
%cellVols: volumes of the cells
%weightDistFactors: Weight distance factors
%distNeighbNods: Distances between cell centers and neigborhood nodes
%solidMask: j true element are over the solid surface

    %Initialize variables*
    u_f=u_f_star;
    v_f=v_f_star;

    %Weast and East Faces*
    for i = 1:ny
        for j=1:nx-1
            
            %Get weight distance factor to east node
            weight_e=weightDistFactors(i,j,3);

            %Get volume from cell P and cell E
            volP=cellVols(i,j);
            volE=cellVols(i,j+1);

            %get distance from node p to node e
            dist_e=distNeighbNods(i,j,3);

            u_f(i,j+1)=u_f_star(i,j+1) + alpha_uv*((volP/aP(i,j))*...
                weight_e + (1-weight_e)*(volE/aP(i,j+1)))*((...
                p_prime(i,j)-p_prime(i,j+1))/dist_e);
        end
    end

    %South and North faces* 
    %interior cells
    for i=1:ny-2
        for j=1:nx

            %Get weight distance factor to south node
            weight_s=weightDistFactors(i,j,4);

            %Get volume from cell P and cell S
            volP=cellVols(i,j);
            volS=cellVols(i+1,j);

            %get distance from node p to node s
            dist_s=distNeighbNods(i,j,4);

            v_f(i+1,j)=v_f_star(i+1,j) + alpha_uv*((volP/aP(i,j))*...
                weight_s + (1-weight_s)*(volS/aP(i+1,j)))*((...
                p_prime(i,j)-p_prime(i+1,j))/dist_s);
            
        end
    end
    
    %south edge*
    i=ny-1;
    for j=1:nx

        %Get weight distance factor to south node
        weight_s=weightDistFactors(i,j,4);

        %Get volume from cell P and cell S
        volP=cellVols(i,j);
        volS=cellVols(i+1,j);

        %get distance from node p to node s
        dist_s=distNeighbNods(i,j,4);

        if solidMask(j) %Solid wall    
            v_f(i+1,j)=v_f_star(i+1,j) + alpha_uv*((volP/aP(i,j))*...
                weight_s + (1-weight_s)*(volS/aP(i+1,j)))*((...
                p_prime(i,j)-p_prime(i+1,j))/dist_s);
        else%Symetric condition (use Apv)
            v_f(i+1,j)=v_f_star(i+1,j) + alpha_uv*((volP/aP(i,j))*...
                weight_s + (1-weight_s)*(volS/aPv(j)))*((...
                p_prime(i,j)-p_prime(i+1,j))/dist_s);
        end
    end
    
end