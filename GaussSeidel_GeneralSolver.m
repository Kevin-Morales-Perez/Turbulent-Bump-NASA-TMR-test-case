    function [phi,rsid_phi,err_phi,iterations_p] = GaussSeidel_GeneralSolver(aP,aW,...
    aN,aE,aS,Su,rsid_phi,epsilon_phi,err_phi,max_iterations_phi,...
    alpha_phi,nx,ny,phi,phi_last)

%GAUSS SEIDEL SOLVER FOR ANY SYSTEM OF LINEAR EQUATIONS WITH 
%IMPLICIT UNDER-RELAXATION APPLIED

%aP Diagonal Coeffitients
%aW,aN,aE,aS: Neigborhood Coeffitients
%Su: Source
%rsid_phi:Array of raw residuals
%epsilon_phi; target error
%err_phi: L2 norm of Raw residuals
%max_iterations_phi: Max inner iterations allowed
%alpha_phi: Under - Relaxation Applied
%nx,ny: size of the mesh
%phi <-- transported varible , used from last it to initialize and 
% accelerate convergence - Changes as iterationd go
%phi_last: phi from last iteration for implicit underelaxation,
%  remain fixed

    
    %Initialice inner iterations counter
    iterations_p=0;
    convergedFlg_phi=false;

    while convergedFlg_phi==false

        iterations_p=iterations_p +1;

        
        %Interior cells
        for i=2:ny-1
            for j=2:nx-1
                phi(i,j)= (alpha_phi/aP(i,j))*(aW(i,j)*phi(i,j-1) + ...
                    aN(i,j)*phi(i-1,j) + aE(i,j)*phi(i,j+1) + ...
                    aS(i,j)*phi(i+1,j) + Su(i,j)) + ...
                    (1-alpha_phi)*phi_last(i,j);
            end
        end

        %EDGES
        %West
        j=1;
        for i=2:ny-1
            phi(i,j)= (alpha_phi/aP(i,j))*(aN(i,j)*phi(i-1,j) + ...
                aE(i,j)*phi(i,j+1) + aS(i,j)*phi(i+1,j) + Su(i,j)) + ...
                (1-alpha_phi)*phi_last(i,j);
        end

        %North
        i=1;
        for j=2:nx-1
            phi(i,j)= (alpha_phi/aP(i,j))*(aW(i,j)*phi(i,j-1) + ...
                aE(i,j)*phi(i,j+1) + aS(i,j)*phi(i+1,j) + Su(i,j)) + ...
                (1-alpha_phi)*phi_last(i,j);

        end

        %East
        j=nx;
        for i=2:ny-1
            phi(i,j)= (alpha_phi/aP(i,j))*(aW(i,j)*phi(i,j-1) + ...
                aN(i,j)*phi(i-1,j) + aS(i,j)*phi(i+1,j) + Su(i,j)) +...
                (1-alpha_phi)*phi_last(i,j);
        end

        %South
        i=ny;
        for j=2:nx-1
            phi(i,j)= (alpha_phi/aP(i,j))*(aW(i,j)*phi(i,j-1) + ...
                aN(i,j)*phi(i-1,j) + aE(i,j)*phi(i,j+1) + Su(i,j)) + ...
                (1-alpha_phi)*phi_last(i,j);

        end

        %Corners

        %West-North
        i=1;
        j=1;

        phi(i,j)= (alpha_phi/aP(i,j))*(aE(i,j)*phi(i,j+1) + ...
            aS(i,j)*phi(i+1,j) + Su(i,j)) + (1-alpha_phi)*phi_last(i,j);

        %East-North
        i=1;
        j=nx;

        phi(i,j)= (alpha_phi/aP(i,j))*(aW(i,j)*phi(i,j-1) + ...
            aS(i,j)*phi(i+1,j) + Su(i,j)) + (1-alpha_phi)*phi_last(i,j);

        %East-South
        i=ny;
        j=nx;

        phi(i,j)= (alpha_phi/aP(i,j))*(aW(i,j)*phi(i,j-1) + ...
            aN(i,j)*phi(i-1,j) + Su(i,j)) + (1-alpha_phi)*phi_last(i,j);

        %West-South
        i=ny;
        j=1;

        phi(i,j)= (alpha_phi/aP(i,j))*(aE(i,j)*phi(i,j+1) + ...
            aN(i,j)*phi(i-1,j) + Su(i,j)) + (1-alpha_phi)*phi_last(i,j);
           
        %Residual computation___________________________________________
        rsid_phi=zeros(ny,nx);% RESET RAW RESIDUALS

        %Interior cells
        for i=2:ny-1
            for j=2:nx-1
                rsid_phi(i,j)= phi(i,j) - ...
                    ((alpha_phi/aP(i,j))*(aW(i,j)*phi(i,j-1) + ...
                    aN(i,j)*phi(i-1,j) + aE(i,j)*phi(i,j+1) + ...
                    aS(i,j)*phi(i+1,j) + Su(i,j)) + ...
                    (1-alpha_phi)*phi_last(i,j));
            end
        end

        %EDGES
        %West
        j=1;
        for i=2:ny-1
            rsid_phi(i,j) = phi(i,j) - ...
                ((alpha_phi/aP(i,j))*(aN(i,j)*phi(i-1,j) + ...
                aE(i,j)*phi(i,j+1) + aS(i,j)*phi(i+1,j) + Su(i,j)) + ...
                (1-alpha_phi)*phi_last(i,j));
        end

        %North
        i=1;
        for j=2:nx-2
            rsid_phi(i,j) = phi(i,j) - ...
                ((alpha_phi/aP(i,j))*(aW(i,j)*phi(i,j-1) + ...
                aE(i,j)*phi(i,j+1) + aS(i,j)*phi(i+1,j) + Su(i,j)) + ...
                (1-alpha_phi)*phi_last(i,j));

        end

        %East
        j=nx;
        for i=2:ny-2
            rsid_phi(i,j)= phi(i,j) - ...
                ((alpha_phi/aP(i,j))*(aW(i,j)*phi(i,j-1) + ...
                aN(i,j)*phi(i-1,j) + aS(i,j)*phi(i+1,j) + Su(i,j)) +...
                (1-alpha_phi)*phi_last(i,j));
        end

        %South
        i=ny;
        for j=2:nx-2
            rsid_phi(i,j) = phi(i,j) - ...
                ((alpha_phi/aP(i,j))*(aW(i,j)*phi(i,j-1) + ...
                aN(i,j)*phi(i-1,j) + aE(i,j)*phi(i,j+1) + Su(i,j)) + ...
                (1-alpha_phi)*phi_last(i,j));

        end

        %Corners

        %West-North
        i=1;
        j=1;

        rsid_phi(i,j) = phi(i,j) - ((alpha_phi/aP(i,j))*(aE(i,j)*phi(i,j+1) + ...
            aS(i,j)*phi(i+1,j) + Su(i,j)) + (1-alpha_phi)*phi_last(i,j));

        %East-North
        i=1;
        j=nx;

        rsid_phi(i,j) = phi(i,j) - ((alpha_phi/aP(i,j))*(aW(i,j)*phi(i,j-1) + ...
            aS(i,j)*phi(i+1,j) + Su(i,j)) + (1-alpha_phi)*phi_last(i,j));

        %East-South
        i=ny;
        j=nx;

        rsid_phi(i,j) = phi(i,j) - ((alpha_phi/aP(i,j))*(aW(i,j)*phi(i,j-1) + ...
            aN(i,j)*phi(i-1,j) + Su(i,j)) + (1-alpha_phi)*phi_last(i,j));

        %West-South
        i=ny;
        j=1;

        rsid_phi(i,j) = phi(i,j) - ((alpha_phi/aP(i,j))*(aE(i,j)*phi(i,j+1) + ...
            aN(i,j)*phi(i-1,j) + Su(i,j)) + (1-alpha_phi)*phi_last(i,j));

        
        %Error from residual L2 norm____________________________________
        err_phi=rms(rsid_phi(:));

        %CHECK CONVERGENCE CRITERIA______________________________________
        if (iterations_p>max_iterations_phi)||(err_phi<epsilon_phi)
            convergedFlg_phi=true;
        end

    end
end