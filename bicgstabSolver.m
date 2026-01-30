function [u,rsid_u,err_u,err_x,total_it] = bicgstabSolver(aP,aW,aN,aE,aS,suX,...
    A_u,suX_vec,restart_u,tol_u,maxit_u,...
    indx_interior,indx_edg_W,indx_edg_N,indx_edg_E,indx_edg_S,...
    indx_cor_WN,indx_cor_EN,indx_cor_ES,indx_cor_WS, ...
    indx_mat,ny,nx,ncells,u_star,...
    alpha_uv)
    %u: solution
    %rsid_u: residuals matrix
    %err_u: l2 norm
    %err_x:output error from GMRES function
    %total_it: Inner *outer iterations in GMRES function 


    u=zeros(ny,nx);
    
    u_star_vec=zeros(ncells,1);


    %INTERIOR CELLS 
    for k=indx_interior
    
        i=indx_mat(k,1);
        j=indx_mat(k,2);
    
        A_u(k,k-1)=-aW(i,j);
        A_u(k,k-nx)=-aN(i,j);
        A_u(k,k+1)=-aE(i,j);
        A_u(k,k+nx)=-aS(i,j);
        A_u(k,k)=aP(i,j)/alpha_uv;
        
    end 
    
    %EDGES
    for k=indx_edg_W %West Edge (4)
    
        i=indx_mat(k,1);
        j=indx_mat(k,2);
    
        A_u(k,k-nx)=-aN(i,j);
        A_u(k,k+1)=-aE(i,j);
        A_u(k,k+nx)=-aS(i,j);
        A_u(k,k)=aP(i,j)/alpha_uv;
    
    end
    for k=indx_edg_N %North edge (2)
    
        i=indx_mat(k,1);
        j=indx_mat(k,2);
    
        A_u(k,k-1)=-aW(i,j);
        A_u(k,k+1)=-aE(i,j);
        A_u(k,k+nx)=-aS(i,j);
        A_u(k,k)=aP(i,j)/alpha_uv;
    
    end
    for k=indx_edg_E %East Edge (5)
    
        i=indx_mat(k,1);
        j=indx_mat(k,2);
    
        A_u(k,k-1)=-aW(i,j);
        A_u(k,k-nx)=-aN(i,j);
        A_u(k,k+nx)=-aS(i,j);
        A_u(k,k)=aP(i,j)/alpha_uv;
    
    end
    for k=indx_edg_S %South Edge (7)
    
        i=indx_mat(k,1);
        j=indx_mat(k,2);
    
        A_u(k,k-1)=-aW(i,j);
        A_u(k,k-nx)=-aN(i,j);
        A_u(k,k+1)=-aE(i,j);
        A_u(k,k)=aP(i,j)/alpha_uv;
    
    end
    
    %CORNERS
    
    k=indx_cor_WN; %West North corner (1)
    
    i=indx_mat(k,1);
    j=indx_mat(k,2);
    
    A_u(k,k+1)=-aE(i,j);
    A_u(k,k+nx)=-aS(i,j);
    A_u(k,k)=aP(i,j)/alpha_uv;
    
    k=indx_cor_EN; %East North corner(3)
    
    i=indx_mat(k,1);
    j=indx_mat(k,2);
    
    A_u(k,k-1)=-aW(i,j);
    A_u(k,k+nx)=-aS(i,j);
    A_u(k,k)=aP(i,j)/alpha_uv;
    
    k=indx_cor_WS; %West South corner (6) 
    
    i=indx_mat(k,1);
    j=indx_mat(k,2);
    
    A_u(k,k-nx)=-aN(i,j);
    A_u(k,k+1)=-aE(i,j);
    A_u(k,k)=aP(i,j)/alpha_uv;
    
    k=indx_cor_ES; %East South corner (8)
    
    i=indx_mat(k,1);
    j=indx_mat(k,2);
    
    A_u(k,k-1)=-aW(i,j);
    A_u(k,k-nx)=-aN(i,j);
    A_u(k,k)=aP(i,j)/alpha_uv;
    
    for k=1:ncells
    
        i=indx_mat(k,1);
        j=indx_mat(k,2);
    
        suX_vec(k,1)=suX(i,j) + ...
            ((1-alpha_uv)/alpha_uv)*A_u(k,k)*u_star(i,j);
        u_star_vec(k,1)=u_star(i,j);
    
    end
    
    %ilu PRECONDITIONER
    [L_u,U_u]=ilu(A_u,struct('type','ilutp','droptol',1e-3,'udiag',1));
    %Solving with GMRS
    [u_vec,flag_n,err_x,iter_n,resvec]=...
    bicgstabl(A_u,suX_vec,restart_u,tol_u,maxit_u,L_u,U_u);

    it_in=iter_n(1);
    it_out=iter_n(2);

    total_it=it_in*it_out;
    
    % Return to matrix form
    for k=1:ncells

        i=indx_mat(k,1);
        j=indx_mat(k,2);
        
        u(i,j)=u_vec(k);

    end

    %L2 NORM
    rsid_u=zeros(ny,nx);% RESET RAW RESIDUALS

    %Interior cells
    for i=2:ny-1
        for j=2:nx-1
            rsid_u(i,j)= u(i,j) - ...
                ((alpha_uv/aP(i,j))*(aW(i,j)*u(i,j-1) + ...
                aN(i,j)*u(i-1,j) + aE(i,j)*u(i,j+1) + ...
                aS(i,j)*u(i+1,j) + suX(i,j)) + ...
                (1-alpha_uv)*u_star(i,j));
        end
    end

    %EDGES
    %West
    j=1;
    for i=2:ny-1
        rsid_u(i,j) = u(i,j) - ...
            ((alpha_uv/aP(i,j))*(aN(i,j)*u(i-1,j) + ...
            aE(i,j)*u(i,j+1) + aS(i,j)*u(i+1,j) + suX(i,j)) + ...
            (1-alpha_uv)*u_star(i,j));
    end

    %North
    i=1;
    for j=2:nx-2
        rsid_u(i,j) = u(i,j) - ...
            ((alpha_uv/aP(i,j))*(aW(i,j)*u(i,j-1) + ...
            aE(i,j)*u(i,j+1) + aS(i,j)*u(i+1,j) + suX(i,j)) + ...
            (1-alpha_uv)*u_star(i,j));

    end

    %East
    j=nx;
    for i=2:ny-2
        rsid_u(i,j)= u(i,j) - ...
            ((alpha_uv/aP(i,j))*(aW(i,j)*u(i,j-1) + ...
            aN(i,j)*u(i-1,j) + aS(i,j)*u(i+1,j) + suX(i,j)) +...
            (1-alpha_uv)*u_star(i,j));
    end

    %South
    i=ny;
    for j=2:nx-2
        rsid_u(i,j) = u(i,j) - ...
            ((alpha_uv/aP(i,j))*(aW(i,j)*u(i,j-1) + ...
            aN(i,j)*u(i-1,j) + aE(i,j)*u(i,j+1) + suX(i,j)) + ...
            (1-alpha_uv)*u_star(i,j));

    end

    %Corners

    %West-North
    i=1;
    j=1;

    rsid_u(i,j) = u(i,j) - ((alpha_uv/aP(i,j))*(aE(i,j)*u(i,j+1) + ...
        aS(i,j)*u(i+1,j) + suX(i,j)) + (1-alpha_uv)*u_star(i,j));

    %East-North
    i=1;
    j=nx;

    rsid_u(i,j) = u(i,j) - ((alpha_uv/aP(i,j))*(aW(i,j)*u(i,j-1) + ...
        aS(i,j)*u(i+1,j) + suX(i,j)) + (1-alpha_uv)*u_star(i,j));

    %East-South
    i=ny;
    j=nx;

    rsid_u(i,j) = u(i,j) - ((alpha_uv/aP(i,j))*(aW(i,j)*u(i,j-1) + ...
        aN(i,j)*u(i-1,j) + suX(i,j)) + (1-alpha_uv)*u_star(i,j));

    %West-South
    i=ny;
    j=1;

    rsid_u(i,j) = u(i,j) - ((alpha_uv/aP(i,j))*(aE(i,j)*u(i,j+1) + ...
        aN(i,j)*u(i-1,j) + suX(i,j)) + (1-alpha_uv)*u_star(i,j));

    
    %Error from residual L2 norm____________________________________
    err_u=rms(rsid_u(:));

    



    for k=1:ncells

        i=indx_mat(k,1);
        j=indx_mat(k,2);
        
        u(i,j)=u_vec(k);

    end
   
end