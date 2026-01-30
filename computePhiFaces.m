function [mu_turbulent_fw,mu_turbulent_fn,mu_turbulent_fe,mu_turbulent_fs]...
    = computePhiFaces(mu_turbulent,weightDistFactors,nx,ny)
    %COMPUTE TRASNPORTED VARIABLE PHI AT FACES USING PHI AT CELL CENTROIDS 
    %AND WEIGHT FACTORS 

    %initialize variables           
    mu_turbulent_fw=zeros(ny,nx);
    mu_turbulent_fn=zeros(ny,nx);
    mu_turbulent_fe=zeros(ny,nx);
    mu_turbulent_fs=zeros(ny,nx);
    
    %West faces
    for i=1:ny
        for j=2:nx
            wf=weightDistFactors(i,j,1);
            mu_turbulent_fw(i,j)=mu_turbulent(i,j)*wf + ...
                mu_turbulent(i,j-1)*(1-wf);
        end
    end
    %First order interpolation for west edge
    mu_turbulent_fw(:,1)=mu_turbulent(:,1);

    %North faces
    for i=2:ny
        for j=1:nx
            wn=weightDistFactors(i,j,2);
            mu_turbulent_fn(i,j)=mu_turbulent(i,j)*wn + ...
                mu_turbulent(i-1,j)*(1-wn);
        end 
    end

    %North edge
    %First order interpolation for north edge
    mu_turbulent_fn(1,:)=mu_turbulent(1,:);


    %East faces
    for i=1:ny
        for j=1:nx-1
            we=weightDistFactors(i,j,3);
            mu_turbulent_fe(i,j)= mu_turbulent(i,j)*we  + ...
                mu_turbulent(i,j+1)*(1-we);
        end 
    end

    %East edge
    %First order interpolation for east edge
    mu_turbulent_fe(:,nx)=mu_turbulent(:,nx);

    
    %South faces
    for i=1:ny-1
        for j=1:nx
            ws=weightDistFactors(i,j,4);
            mu_turbulent_fs(i,j)=mu_turbulent(i,j)*ws + ...
                mu_turbulent(i+1,j)*(1-ws);
        end
    end

    %South edge
    %First order interpolation for south edge
    mu_turbulent_fs(ny,:)=mu_turbulent(ny,:);



    
