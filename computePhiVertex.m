function [phiVert] = computePhiVertex(phi,phiGradients,...
    VecCentVertx,nx,ny)
%Computation of transported variable phi
%  at vertex of each cell using gradients
%Can be used for Temperature , Velocity and less common Pressure 
%used for nu~ Spallart Allmaras transported variable

    %Initialze variable

    phiVert=zeros(ny,nx,1,4);

    for i=1:ny
        for j=1:nx

            %get gradient from cell i , j
            %grad_phi=reshape(phiGradients(i,j,:,:),[1,2]);

            dphi_x=phiGradients(i,j,1,1);
            dphi_y=phiGradients(i,j,1,2);

            grad_phi=[dphi_x,dphi_y];

            %get temperature from cell i, j

            phi_cell=phi(i,j);

            %get vectors from all vertex of cell i , j

            vecvert=reshape(VecCentVertx(i,j,:,:),[4,2]);

            %split vectors

            vec_wn=vecvert(1,:);
            vec_en=vecvert(2,:);
            vec_es=vecvert(3,:);
            vec_ws=vecvert(4,:);

            %compute temperatures at vertex

            phi_wn=phi_cell + dot(grad_phi,vec_wn);

            phi_en=phi_cell + dot(grad_phi,vec_en);

            phi_es=phi_cell + dot(grad_phi,vec_es);

            phi_ws=phi_cell + dot(grad_phi,vec_ws);

            %Store values 

            phiVert(i,j,:,:)=[phi_wn,phi_en,phi_es,phi_ws];

        end
    end

end
