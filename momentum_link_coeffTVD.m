function [aW,aE,aN,aS,aP,aPv,suX,suY] = momentum_link_coeffTVD(rho,mu,...
    mu_turbulent_fw,mu_turbulent_fn,mu_turbulent_fe,mu_turbulent_fs,...
    lgtFaces,cellVols,dW,dN,dE,dS,dW_c,dN_c,dE_c,dS_c,u0,u_f,v_f,...
    u_cor,v_cor,grad_p,nx,ny,solidMask,u,v,grad_u,grad_v,VecCentNbNods)
%COEFFITIENTS FOR MOMENTUM SYSTEM EQUATIONS 

%INITIALIZE VARIABLES
aW=zeros(ny,nx);
aN=zeros(ny,nx);
aE=zeros(ny,nx);
aS=zeros(ny,nx);
aP=zeros(ny,nx);
aPv=zeros(1,nx);%Diagonal Coeffitienst for simetric boundary condition

%Sources 
suX=zeros(ny,nx);
suY=zeros(ny,nx);

%nx,ny: Mesh sizes
%solidMask: Mask for cells that are marked as above solid surface
%rho:density
%lgtFaces: Area of each face(lenghts in this case because is 2D)
%u_f: velocity normal to faces weast and east in the face centers
%v_f:velocity normal to faces north and south in the face centers
%u0: inlet velocity
%dW,dE,dN,dS:direct diffusion coeffitients
%dW_c,dN_c,dE_c,dS_c:Cross diffusion coeffitients
%aP,aPv:Central coeffitients
%suX,suY: Sources
%cellVols:Volume of each cell
%u_corners,v_corners: velocities at corners
%grad_p:Gradient of pressure 
%grad_tau_mn: Gradient of turbulent Stresses
%mu: dynamic viscosity
%mu_turbulent: Eddy viscosity


%Coeffitients for momentum equations 
%   Convection + Diffusion + Sources Coeffitiens for the momentum equations
%in order to have them in this format considering trasported
%  variable phi:
% Ap*phi_p = Aw*phi_w + An*phi_n + Ae*phi_e + As*phi_s + Su 

        %BOUNDARY  CONDITIONS

%       BOUNDARY CONDITIONS FOR FLOW OVER A BUMP
%           NEUMMAN  (DU/DY=0 ,DV/DY=0)
%                                                      O N  
%I                                                     U E   
%N                                                     T U  DU/DY=0
%L                                                     L M
%E                                                     E M  DV/DY=0
%T                                                     T A 
% ___ SYMETRIC ___|___   NO SLIP ___|___ SYMETRIC  ____  N

%nx_upstr,nx_fp,nx_dwnstr

    %_____________________________________________________________________
    %multipling each difussive term by dynamics viscosity plus eddy visc. 

    mu_total_fw=mu+ mu_turbulent_fw;
    mu_total_fn=mu+ mu_turbulent_fn;
    mu_total_fe=mu+ mu_turbulent_fe;
    mu_total_fs=mu+ mu_turbulent_fs;
    
    dW=dW.*(mu_total_fw);
    dN=dN.*(mu_total_fn);
    dE=dE.*(mu_total_fe);
    dS=dS.*(mu_total_fs);
    dW_c=dW_c.*(mu_total_fw);
    dN_c=dN_c.*(mu_total_fn);    
    dE_c=dE_c.*(mu_total_fe);
    dS_c=dS_c.*(mu_total_fs);

    %_____________________________________________________________________


%1.- MOMENTUM LINK COEFFITIENTS
    
    %1.1- Interior Cells
    for i=2:ny-1
        for j=2:nx-1
            
            %get face areas
            lgt_fw=lgtFaces(i,j,1);
            lgt_fn=lgtFaces(i,j,2); 
            lgt_fe=lgtFaces(i,j,3); 
            lgt_fs=lgtFaces(i,j,4);
  
            %Mass fluxes at faces 
            fW=-rho*lgt_fw*u_f(i,j);
            fN=-rho*lgt_fn*v_f(i,j); 
            fE=rho*lgt_fe*u_f(i,j+1); 
            fS=rho*lgt_fs*v_f(i+1,j);

            %coeffitients using upwind squeme
            aW(i,j)=dW(i,j) + max(0,-fW);
            aE(i,j)=dE(i,j) + max(0,-fE);
            aN(i,j)=dN(i,j) + max(0,-fN);
            aS(i,j)=dS(i,j) + max(0,-fS);
            aP(i,j)=aW(i,j) + aE(i,j) + aN(i,j) + aS(i,j) + fW + fE + fN + fS;

            %explicit computation of cross - difussion components*
            
            %______________________U__________________________
            
            %Get velocity vertexes values
            u_wn=u_cor(i,j,1,1);
            u_en=u_cor(i,j,1,2);
            u_es=u_cor(i,j,1,3);
            u_ws=u_cor(i,j,1,4);

            %compute cross diffusion terms for each face 
            Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
            Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
            Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
            Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
            
            %Sum all contributions
            Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

            %______________________V__________________________

            %Get velocity vertexes values
            v_wn=v_cor(i,j,1,1);
            v_en=v_cor(i,j,1,2);
            v_es=v_cor(i,j,1,3);
            v_ws=v_cor(i,j,1,4);

            %compute cross diffusion terms for each face 
            Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
            Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
            Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
            Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
            
            %Sum all contributions
            Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;
                              
            %____________Get pressure Gradients* __________________
            
            %get the gradient of pressure at the center of the cell
            
            %Grad in x
            grad_px=grad_p(i,j,1,1);
            %Grad in y
            grad_py=grad_p(i,j,1,2);

            %__________________TVD Corrections___________________________ 

            %Vector with mass fluxes 
            f_C=[fW,fN,fE,fS];

            %Vector with nodal values
            u_vec=[u(i,j-1),u(i-1,j),u(i,j+1),u(i+1,j),u(i,j)];
            v_vec=[v(i,j-1),v(i-1,j),v(i,j+1),v(i+1,j),v(i,j)];

            %vector with Gradients for U
            grad_u_mat=[grad_u(i,j-1,1,1),grad_u(i,j-1,1,2);...
                grad_u(i-1,j,1,1),grad_u(i-1,j,1,2);...
                grad_u(i,j+1,1,1),grad_u(i,j+1,1,2);...
                grad_u(i+1,j,1,1),grad_u(i+1,j,1,2);...
                grad_u(i,j,1,1),grad_u(i,j,1,2)];

            %vector with gradients for V
            grad_v_mat=[grad_v(i,j-1,1,1),grad_v(i,j-1,1,2);...
                grad_v(i-1,j,1,1),grad_v(i-1,j,1,2);...
                grad_v(i,j+1,1,1),grad_v(i,j+1,1,2);...
                grad_v(i+1,j,1,1),grad_v(i+1,j,1,2);...
                grad_v(i,j,1,1),grad_v(i,j,1,2)];

            %Vector with vectors from P to Neigborhood nodes
            r_mat=[VecCentNbNods(i,j,1,1),VecCentNbNods(i,j,1,2);...
                VecCentNbNods(i,j,2,1),VecCentNbNods(i,j,2,2);...
                VecCentNbNods(i,j,3,1),VecCentNbNods(i,j,3,2);...
                VecCentNbNods(i,j,4,1),VecCentNbNods(i,j,4,2)];


            [s_dc_total_x,s_dc_total_y] = tvdDefCorr(f_C,u_vec,v_vec,...
                grad_u_mat,grad_v_mat,r_mat);


            %_____ Sources (Pressure derivative + cross difussion + tvd)* ______

            suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + s_dc_total_x;
            suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + s_dc_total_y;
    
        end
    end
    
    %1.2- Walls
    
    %1.2.1- Left wall - West (INLET)*
    
    j=1;
    
    for i =2:ny-1

        %get face areas
        lgt_fw=lgtFaces(i,j,1);
        lgt_fn=lgtFaces(i,j,2);
        lgt_fe=lgtFaces(i,j,3);
        lgt_fs=lgtFaces(i,j,4);

        %Mass fluxes at faces 
        fW=-rho*lgt_fw*u_f(i,j);
        fN=-rho*lgt_fn*v_f(i,j); 
        fE=rho*lgt_fe*u_f(i,j+1);         
        fS=rho*lgt_fs*v_f(i+1,j);   
            
        %coeffitients using upwind squeme
        aW(i,j)=0;%No adyacent w node  
        aE(i,j)=dE(i,j) + max(0,-fE);
        aN(i,j)=dN(i,j) + max(0,-fN);
        aS(i,j)=dS(i,j) + max(0,-fS);
        aP(i,j)=dW(i,j) + aE(i,j) + aN(i,j) + aS(i,j) +fN +fE +fS;

        %____________ CROSS DIFFUSION TERMS______________________
        %______________________U__________________________
            
        %Get velocity vertexes values
        

        %Split velocities
        u_wn=u_cor(i,j,1,1);
        u_en=u_cor(i,j,1,2);
        u_es=u_cor(i,j,1,3);
        u_ws=u_cor(i,j,1,4);

        %compute cross diffusion terms for each face 
        Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
        Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
        Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
        Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
        
        %Sum all contributions
        Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

            %______________________V__________________________

        %Get velocity vertexes values
        
        %Split velocities
        v_wn=v_cor(i,j,1,1);
        v_en=v_cor(i,j,1,2);
        v_es=v_cor(i,j,1,3);
        v_ws=v_cor(i,j,1,4);

        %compute cross diffusion terms for each face 
        Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
        Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
        Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
        Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
        
        %Sum all contributions
        Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;

        %____________Get pressure Gradients__________________
            
        %get the gradient of pressure at the center of the cell
        
        
        %Grad in x
        grad_px=grad_p(i,j,1,1);
        %Grad in y
        grad_py=grad_p(i,j,1,2);

        %__________________TVD Corrections___________________________ 

        %Vector with mass fluxes 
        f_C=[fN,fE,fS];

        %Vector with nodal values
        u_vec=[u(i-1,j),u(i,j+1),u(i+1,j),u(i,j)];
        v_vec=[v(i-1,j),v(i,j+1),v(i+1,j),v(i,j)];

        %vector with Gradients for U
        grad_u_mat=[grad_u(i-1,j,1,1),grad_u(i-1,j,1,2);...
            grad_u(i,j+1,1,1),grad_u(i,j+1,1,2);...
            grad_u(i+1,j,1,1),grad_u(i+1,j,1,2);...
            grad_u(i,j,1,1),grad_u(i,j,1,2)];

        %vector with gradients for V
        grad_v_mat=[grad_v(i-1,j,1,1),grad_v(i-1,j,1,2);...
            grad_v(i,j+1,1,1),grad_v(i,j+1,1,2);...
            grad_v(i+1,j,1,1),grad_v(i+1,j,1,2);...
            grad_v(i,j,1,1),grad_v(i,j,1,2)];

        %Vector with vectors from P to Neigborhood nodes
        r_mat=[VecCentNbNods(i,j,2,1),VecCentNbNods(i,j,2,2);...
            VecCentNbNods(i,j,3,1),VecCentNbNods(i,j,3,2);...
            VecCentNbNods(i,j,4,1),VecCentNbNods(i,j,4,2)];


        [s_dc_total_x,s_dc_total_y] = tvdDefCorr(f_C,u_vec,v_vec,...
            grad_u_mat,grad_v_mat,r_mat);

        %Sources (Pressure derivative + TVD)
        suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + u0*(dW(i,j) -fW) + ...
            s_dc_total_x;  
        suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + s_dc_total_y;
    
    end
    
    %1.2.2- Top wall - North  (NEUMMAN)*
    
    i=1;

    for j=2:nx-1

        %get face areas
        lgt_fw=lgtFaces(i,j,1);
        lgt_fn=lgtFaces(i,j,2);%Neumman condition
        lgt_fe=lgtFaces(i,j,3);
        lgt_fs=lgtFaces(i,j,4);
    
        %Mass fluxes
        fW=-rho*lgt_fw*u_f(i,j);
        fN=-rho*lgt_fn*v_f(i,j); 
        fE=rho*lgt_fe*u_f(i,j+1);        
        fS=rho*lgt_fs*v_f(i+1,j);   
            
        %coeffitients using upwind scheme
        aW(i,j)=dW(i,j) + max(0,-fW);
        aE(i,j)=dE(i,j) + max(0,-fE);
        aN(i,j)=0;%No adyacent node N
        aS(i,j)=dS(i,j) + max(0,-fS);
        aP(i,j)=aW(i,j) + aE(i,j)  + aS(i,j) + fW + fN + fE + fS;


        %explicit computation of cross - difussion components            
        %______________________U__________________________
        
        %Get velocity vertexes values
        

        %Split velocities
        u_wn=u_cor(i,j,1,1);
        u_en=u_cor(i,j,1,2);
        u_es=u_cor(i,j,1,3);
        u_ws=u_cor(i,j,1,4);

        %compute cross diffusion terms for each face 
        Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
        Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
        Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
        Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
        
        %Sum all contributions
        Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;
        
        %______________________v__________________________

        %Get velocity vertexes values
        
        %Split velocities
        v_wn=v_cor(i,j,1,1);
        v_en=v_cor(i,j,1,2);
        v_es=v_cor(i,j,1,3);
        v_ws=v_cor(i,j,1,4);

        %compute cross diffusion terms for each face 
        Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
        Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
        Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
        Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
        
        %Sum all contributions
        Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;
        
        %____________Compute pressure Gradients__________________

        %get the gradient of pressure at the center of the cell

        %Grad in x
        grad_px=grad_p(i,j,1,1);
        %Grad in y
        grad_py=grad_p(i,j,1,2);

        %__________________TVD Corrections___________________________ 

        %Vector with mass fluxes 
        f_C=[fW,fE,fS];

        %Vector with nodal values
        u_vec=[u(i,j-1),u(i,j+1),u(i+1,j),u(i,j)];
        v_vec=[v(i,j-1),v(i,j+1),v(i+1,j),v(i,j)];

        %vector with Gradients for U
        grad_u_mat=[grad_u(i,j-1,1,1),grad_u(i,j-1,1,2);...
            grad_u(i,j+1,1,1),grad_u(i,j+1,1,2);...
            grad_u(i+1,j,1,1),grad_u(i+1,j,1,2);...
            grad_u(i,j,1,1),grad_u(i,j,1,2)];

        %vector with gradients for V
        grad_v_mat=[grad_v(i,j-1,1,1),grad_v(i,j-1,1,2);...
            grad_v(i,j+1,1,1),grad_v(i,j+1,1,2);...
            grad_v(i+1,j,1,1),grad_v(i+1,j,1,2);...
            grad_v(i,j,1,1),grad_v(i,j,1,2)];

        %Vector with vectors from P to Neigborhood nodes
        r_mat=[VecCentNbNods(i,j,1,1),VecCentNbNods(i,j,1,2);...
            VecCentNbNods(i,j,3,1),VecCentNbNods(i,j,3,2);...
            VecCentNbNods(i,j,4,1),VecCentNbNods(i,j,4,2)];


        [s_dc_total_x,s_dc_total_y] = tvdDefCorr(f_C,u_vec,v_vec,...
            grad_u_mat,grad_v_mat,r_mat);

        %Sources (Pressure derivative + boundary condition)
        suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + s_dc_total_x; 
        suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + s_dc_total_y;
    
    end
    
    %1.2.3- Right Wall - East  (Outlet)*
    
    j=nx;
    for i=2:ny-1
       
        %get face areas
        lgt_fw=lgtFaces(i,j,1);
        lgt_fn=lgtFaces(i,j,2);
        lgt_fe=lgtFaces(i,j,3);
        lgt_fs=lgtFaces(i,j,4);

        %Mass fluxes at faces
        fW=-rho*lgt_fw*u_f(i,j);
        fN=-rho*lgt_fn*v_f(i,j); 
        fE=rho*lgt_fe*u_f(i,j+1); %Neumman condition         
        fS=rho*lgt_fs*v_f(i+1,j);   
        
        %coeffitients using upwind squeme
        aW(i,j)=dW(i,j) + max(0,-fW);
        aE(i,j)=0;
        aN(i,j)=dN(i,j) + max(0,-fN);
        aS(i,j)=dS(i,j) + max(0,-fS);
        aP(i,j)=aW(i,j) + aN(i,j) + aS(i,j)  +fW +fN +fS +fE;

        %____________  CROSS DIFFUSION TERMS______________________

        %______________________U__________________________
            
        %Get velocity vertexes values
        

        %Split velocities
        u_wn=u_cor(i,j,1,1);
        u_en=u_cor(i,j,1,2);
        u_es=u_cor(i,j,1,3);
        u_ws=u_cor(i,j,1,4);

        %compute cross diffusion terms for each face 
        Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
        Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
        Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
        Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
        
        %Sum all contributions
        Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

        %______________________V__________________________

        %Get velocity vertexes values
        

        %Split velocities
        v_wn=v_cor(i,j,1,1);
        v_en=v_cor(i,j,1,2);
        v_es=v_cor(i,j,1,3);
        v_ws=v_cor(i,j,1,4);

        %compute cross diffusion terms for each face 
        Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
        Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
        Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
        Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
        
        %Sum all contributions
        Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;

        %____________Get pressure Gradients__________________
        
        %get the gradient of pressure at the center of the cell
        
        %Grad in x
        grad_px=grad_p(i,j,1,1);
        %Grad in y
        grad_py=grad_p(i,j,1,2);

        %__________________TVD Corrections___________________________ 

        %Vector with mass fluxes 
        f_C=[fW,fN,fS];

        %Vector with nodal values
        u_vec=[u(i,j-1),u(i-1,j),u(i+1,j),u(i,j)];
        v_vec=[v(i,j-1),v(i-1,j),v(i+1,j),v(i,j)];

        %vector with Gradients for U
        grad_u_mat=[grad_u(i,j-1,1,1),grad_u(i,j-1,1,2);...
            grad_u(i-1,j,1,1),grad_u(i-1,j,1,2);...
            grad_u(i+1,j,1,1),grad_u(i+1,j,1,2);...
            grad_u(i,j,1,1),grad_u(i,j,1,2)];

        %vector with gradients for V
        grad_v_mat=[grad_v(i,j-1,1,1),grad_v(i,j-1,1,2);...
            grad_v(i-1,j,1,1),grad_v(i-1,j,1,2);...
            grad_v(i+1,j,1,1),grad_v(i+1,j,1,2);...
            grad_v(i,j,1,1),grad_v(i,j,1,2)];

        %Vector with vectors from P to Neigborhood nodes
        r_mat=[VecCentNbNods(i,j,1,1),VecCentNbNods(i,j,1,2);...
            VecCentNbNods(i,j,2,1),VecCentNbNods(i,j,2,2);...
            VecCentNbNods(i,j,4,1),VecCentNbNods(i,j,4,2)];


        [s_dc_total_x,s_dc_total_y] = tvdDefCorr(f_C,u_vec,v_vec,...
            grad_u_mat,grad_v_mat,r_mat);


        %_____ Sources (Pressure derivative + cross difussion) ______

        suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + s_dc_total_x;
        suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + s_dc_total_y;
    
    end
    
    %1.2.4- Bottom wall -- South (Symetric in the free-stream and no slip 
    % on the plate)*
    
    i=ny;
    
    for j=2:nx-1

        %get face areas
        lgt_fw=lgtFaces(i,j,1);
        lgt_fn=lgtFaces(i,j,2); 
        lgt_fe=lgtFaces(i,j,3); 
        %lgt_fs=lgtFaces(i,j,4);

        %Mass fluxes at faces 
        fW=-rho*lgt_fw*u_f(i,j);
        fN=-rho*lgt_fn*v_f(i,j);
        fE=rho*lgt_fe*u_f(i,j+1);  
        %fS=rho*lgt_fs*v_f(i+1,j); No mass flux at south face

        %coeffitients using upwind squeme
        aW(i,j)=dW(i,j) + max(0,-fW);
        aE(i,j)=dE(i,j) + max(0,-fE);
        aN(i,j)=dN(i,j) + max(0,-fN);
        %aS(i,j)=0;

        if solidMask(j) %Solid part (Dirichelet , fixed 0 for u and v)
            %central Coeffitient for u and v are the same (over the plate)
            aP(i,j)=aW(i,j) + aE(i,j) + aN(i,j) + dS(i,j) +fW +fN +fE;
          
        else %free part %(Neumman fo u , dirichlet fixed 0 for v)

            %central coeffitient for momentum equation for v will be
            %different due to Symetric BC (free- stream zones)
            aP(i,j)=aW(i,j) + aE(i,j) + aN(i,j) +fW +fN +fE;
            aPv(j)=aW(i,j) + aE(i,j) + aN(i,j)+ dS(i,j) +fW +fN +fE;
 
        end

        %explicit computation of cross - difussion components
        %______________________U__________________________
        
        %Get velocity vertexes values
        
        %Split velocities
        u_wn=u_cor(i,j,1,1);
        u_en=u_cor(i,j,1,2);
        u_es=u_cor(i,j,1,3);
        u_ws=u_cor(i,j,1,4);

        %compute cross diffusion terms for each face 
        Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
        Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
        Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
        Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
        
        %Sum all contributions
        Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

        %______________________V__________________________

        %Get velocity vertexes values
        
        %Split velocities
        v_wn=v_cor(i,j,1,1);
        v_en=v_cor(i,j,1,2);
        v_es=v_cor(i,j,1,3);
        v_ws=v_cor(i,j,1,4);

        %compute cross diffusion terms for each face 
        Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
        Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
        Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
        Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
        
        %Sum all contributions
        Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;
        
        %____________Get pressure Gradients__________________
        
        %get the gradient of pressure at the center of the cell
          
        %Grad in x
        grad_px=grad_p(i,j,1,1);
        %Grad in y
        grad_py=grad_p(i,j,1,2);

        %__________________TVD Corrections___________________________ 

        %Vector with mass fluxes 
        f_C=[fW,fN,fE];

        %Vector with nodal values
        u_vec=[u(i,j-1),u(i-1,j),u(i,j+1),u(i,j)];
        v_vec=[v(i,j-1),v(i-1,j),v(i,j+1),v(i,j)];

        %vector with Gradients for U
        grad_u_mat=[grad_u(i,j-1,1,1),grad_u(i,j-1,1,2);...
            grad_u(i-1,j,1,1),grad_u(i-1,j,1,2);...
            grad_u(i,j+1,1,1),grad_u(i,j+1,1,2);...
            grad_u(i,j,1,1),grad_u(i,j,1,2)];

        %vector with gradients for V
        grad_v_mat=[grad_v(i,j-1,1,1),grad_v(i,j-1,1,2);...
            grad_v(i-1,j,1,1),grad_v(i-1,j,1,2);...
            grad_v(i,j+1,1,1),grad_v(i,j+1,1,2);...
            grad_v(i,j,1,1),grad_v(i,j,1,2)];

        %Vector with vectors from P to Neigborhood nodes
        r_mat=[VecCentNbNods(i,j,1,1),VecCentNbNods(i,j,1,2);...
            VecCentNbNods(i,j,2,1),VecCentNbNods(i,j,2,2);...
            VecCentNbNods(i,j,3,1),VecCentNbNods(i,j,3,2)];

        [s_dc_total_x,s_dc_total_y] = tvdDefCorr(f_C,u_vec,v_vec,...
            grad_u_mat,grad_v_mat,r_mat);

        %Sources (Pressure gradient & TVD corrections)

        suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + s_dc_total_x;
        suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + s_dc_total_y;

    end
    
    %1.3- Corners
    
    %1.3.1- *North-West *************************
    i=1;
    j=1;

    %get face areas
    lgt_fw=lgtFaces(i,j,1);
    lgt_fn=lgtFaces(i,j,2);  
    lgt_fe=lgtFaces(i,j,3); 
    lgt_fs=lgtFaces(i,j,4);
    
    %Mass fluxes at faces 
    fW=-rho*lgt_fw*u_f(i,j);
    fN=-rho*lgt_fn*v_f(i,j); 
    fE=rho*lgt_fe*u_f(i,j+1);    
    fS=rho*lgt_fs*v_f(i+1,j);
    
    %coeffitients using upwind squeme
    aW(i,j)=0;% No adyacent node W
    aE(i,j)=dE(i,j) + max(0,-fE);
    aN(i,j)=0;%No adyacent node N
    aS(i,j)=dS(i,j) + max(0,-fS);
    aP(i,j)=dW(i,j) + aE(i,j) + aS(i,j) + fN + fE +fS;
    
    %____________ CROSS DIFFUSION TERMS______________________
    
    %______________________U__________________________
            
    %Get velocity vertexes values
    

    %Split velocities
    u_wn=u_cor(i,j,1,1);
    u_en=u_cor(i,j,1,2);
    u_es=u_cor(i,j,1,3);
    u_ws=u_cor(i,j,1,4);

    %compute cross diffusion terms for each face 
    Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
    Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
    Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
    Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
    
    %Sum all contributions
    Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

    %______________________V__________________________

    %Get velocity vertexes values
    

    %Split velocities
    v_wn=v_cor(i,j,1,1);
    v_en=v_cor(i,j,1,2);
    v_es=v_cor(i,j,1,3);
    v_ws=v_cor(i,j,1,4);

    %compute cross diffusion terms for each face 
    Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
    Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
    Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
    Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
    
    %Sum all contributions
    Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;

    %____________Get pressure Gradients__________________

    %get the gradient of pressure at the center of the cell
    
    %Grad in x
    grad_px=grad_p(i,j,1,1);
    %Grad in y
    grad_py=grad_p(i,j,1,2);

    %Vector with mass fluxes 
    f_C=[fE,fS];

    %Vector with nodal values
    u_vec=[u(i,j+1),u(i+1,j),u(i,j)];
    v_vec=[v(i,j+1),v(i+1,j),v(i,j)];

    %vector with Gradients for U
    grad_u_mat=[grad_u(i,j+1,1,1),grad_u(i,j+1,1,2);...
        grad_u(i+1,j,1,1),grad_u(i+1,j,1,2);...
        grad_u(i,j,1,1),grad_u(i,j,1,2)];

    %vector with gradients for V
    grad_v_mat=[grad_v(i,j+1,1,1),grad_v(i,j+1,1,2);...
        grad_v(i+1,j,1,1),grad_v(i+1,j,1,2);...
        grad_v(i,j,1,1),grad_v(i,j,1,2)];

    %Vector with vectors from P to Neigborhood nodes
    r_mat=[VecCentNbNods(i,j,3,1),VecCentNbNods(i,j,3,2);...
        VecCentNbNods(i,j,4,1),VecCentNbNods(i,j,4,2)];


    [s_dc_total_x,s_dc_total_y] = tvdDefCorr(f_C,u_vec,v_vec,...
        grad_u_mat,grad_v_mat,r_mat);

    %Sources (Pressure dErivative + boundary condition)
    suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + u0*(dW(i,j) -fW) + s_dc_total_x;
    suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + s_dc_total_y;
    
    %1.3.2- *North-East ************************
    i=1;
    j=nx;
    
    %get face areas
    lgt_fw=lgtFaces(i,j,1); 
    lgt_fn=lgtFaces(i,j,2);
    lgt_fe=lgtFaces(i,j,3);  
    lgt_fs=lgtFaces(i,j,4);

    %Mass fluxes at faces 
    fW=-rho*lgt_fw*u_f(i,j);
    fN=-rho*lgt_fn*v_f(i,j);
    fE=rho*lgt_fe*u_f(i,j+1);     
    fS=rho*lgt_fs*v_f(i+1,j); 
    
    %coeffitients using upwind squeme
    aW(i,j)=dW(i,j) + max(0,-fW);
    aE(i,j)=0;
    aN(i,j)=0;
    aS(i,j)=dS(i,j) + max(0,-fS);
    aP(i,j)=aW(i,j) + aS(i,j) +fW + fN + fS +fE;

    %____________ CROSS DIFFUSION TERMS______________________

    %______________________U__________________________
            
    %Get velocity vertexes values
    

    %Split velocities
    u_wn=u_cor(i,j,1,1);
    u_en=u_cor(i,j,1,2);
    u_es=u_cor(i,j,1,3);
    u_ws=u_cor(i,j,1,4);

    %compute cross diffusion terms for each face 
    Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
    Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
    Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
    Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
    
    %Sum all contributions
    Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

    %______________________V__________________________

    %Get velocity vertexes values
    

    %Split velocities
    v_wn=v_cor(i,j,1,1);
    v_en=v_cor(i,j,1,2);
    v_es=v_cor(i,j,1,3);
    v_ws=v_cor(i,j,1,4);

    %compute cross diffusion terms for each face 
    Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
    Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
    Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
    Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
    
    %Sum all contributions
    Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;

    %____________Get pressure Gradients__________________

    %get the gradient of pressure at the center of the cell
    
    %Grad in x
    grad_px=grad_p(i,j,1,1);
    %Grad in y
    grad_py=grad_p(i,j,1,2);

    %___________________________TVD  Corrections_________________________
    %Vector with mass fluxes 
    f_C=[fW,fS];

    %Vector with nodal values
    u_vec=[u(i,j-1),u(i+1,j),u(i,j)];
    v_vec=[v(i,j-1),v(i+1,j),v(i,j)];

    %vector with Gradients for U
    grad_u_mat=[grad_u(i,j-1,1,1),grad_u(i,j-1,1,2);...
        grad_u(i+1,j,1,1),grad_u(i+1,j,1,2);...
        grad_u(i,j,1,1),grad_u(i,j,1,2)];

    %vector with gradients for V
    grad_v_mat=[grad_v(i,j-1,1,1),grad_v(i,j-1,1,2);...
        grad_v(i+1,j,1,1),grad_v(i+1,j,1,2);...
        grad_v(i,j,1,1),grad_v(i,j,1,2)];

    %Vector with vectors from P to Neigborhood nodes
    r_mat=[VecCentNbNods(i,j,1,1),VecCentNbNods(i,j,1,2);...
        VecCentNbNods(i,j,4,1),VecCentNbNods(i,j,4,2)];


        [s_dc_total_x,s_dc_total_y] = tvdDefCorr(f_C,u_vec,v_vec,...
            grad_u_mat,grad_v_mat,r_mat);

    %Sources (Pressure derivative + boundary condition)
    suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + s_dc_total_x;
    suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + s_dc_total_y;
    
    %1.3.3- *South-East (Symetric at South and outlet at East)***********
    i=ny;
    j=nx;

    %get face areas
    lgt_fw=lgtFaces(i,j,1);
    lgt_fn=lgtFaces(i,j,2);
    lgt_fe=lgtFaces(i,j,3);  
    %lgt_fs=lgtFaces(i,j,4);
    
    %Mass fluxes at faces 
    fW=-rho*lgt_fw*u_f(i,j);
    fN=-rho*lgt_fn*v_f(i,j); 
    fE=rho*lgt_fe*u_f(i,j+1);    
    %fS=rho*lgt_fs*v_f(i+1,j); No mass flux at south face
    
    %coeffitients using upwind scheme
    aW(i,j)=dW(i,j) + max(0,-fW);
    aE(i,j)=0;%No east adyacent node
    aN(i,j)=dN(i,j) + max(0,-fN);
    aS(i,j)=0;%No south adyacent node
    aP(i,j)=aW(i,j) + aN(i,j) + fW + fN + fE;
    aPv(j)= aW(i,j) + aN(i,j) + dS(i,j) +fW +fN +fE;

    %____________ CROSS DIFFUSION TERMS______________________

    %______________________U__________________________
            
    %Get velocity vertexes values
    
    %Split velocities
    u_wn=u_cor(i,j,1,1);
    u_en=u_cor(i,j,1,2);
    u_es=u_cor(i,j,1,3);
    u_ws=u_cor(i,j,1,4);

    %compute cross diffusion terms for each face 
    Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
    Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
    Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
    Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
    
    %Sum all contributions
    Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

    %______________________V__________________________

    %Get velocity vertexes values

    v_wn=v_cor(i,j,1,1);
    v_en=v_cor(i,j,1,2);
    v_es=v_cor(i,j,1,3);
    v_ws=v_cor(i,j,1,4);

    %compute cross diffusion terms for each face 
    Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
    Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
    Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
    Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
    
    %Sum all contributions
    Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;

    %____________Get pressure Gradients__________________

    %get the gradient of pressure at the center of the cell

    %Grad in x
    grad_px=grad_p(i,j,1,1);
    %Grad in y
    grad_py=grad_p(i,j,1,2);

    %__________________TVD Corrections___________________________ 

    %Vector with mass fluxes 
    f_C=[fW,fN];

    %Vector with nodal values
    u_vec=[u(i,j-1),u(i-1,j),u(i,j)];
    v_vec=[v(i,j-1),v(i-1,j),v(i,j)];

    %vector with Gradients for U
    grad_u_mat=[grad_u(i,j-1,1,1),grad_u(i,j-1,1,2);...
        grad_u(i-1,j,1,1),grad_u(i-1,j,1,2);...
        grad_u(i,j,1,1),grad_u(i,j,1,2)];

    %vector with gradients for V
    grad_v_mat=[grad_v(i,j-1,1,1),grad_v(i,j-1,1,2);...
        grad_v(i-1,j,1,1),grad_v(i-1,j,1,2);...
        grad_v(i,j,1,1),grad_v(i,j,1,2)];

    %Vector with vectors from P to Neigborhood nodes
    r_mat=[VecCentNbNods(i,j,1,1),VecCentNbNods(i,j,1,2);...
        VecCentNbNods(i,j,2,1),VecCentNbNods(i,j,2,2)];


    [s_dc_total_x,s_dc_total_y] = tvdDefCorr(f_C,u_vec,v_vec,...
        grad_u_mat,grad_v_mat,r_mat);

    %Sources (Pressure derivative)
    suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + s_dc_total_x;
    suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + s_dc_total_y;
    
    %1.3.4- *South-West (Symetric at South and inlet at West) *************
    i=ny;
    j=1;

    %get face areas
    lgt_fw=lgtFaces(i,j,1);
    lgt_fn=lgtFaces(i,j,2);
    lgt_fe=lgtFaces(i,j,3); 
    %lgt_fs=lgtFaces(i,j,4);
    
    %Mass fluxes at faces 
    fW=-rho*lgt_fw*u_f(i,j);
    fN=-rho*lgt_fn*v_f(i,j); 
    fE=rho*lgt_fe*u_f(i,j+1); 
    %fS=rho*lgt_fs*v_f(i+1,j);
    
    %coeffitients using upwind squeme
    aW(i,j)=0;
    aE(i,j)=dE(i,j) + max(0,-fE);
    aN(i,j)=dN(i,j) + max(0,-fN);
    aS(i,j)=0;
    aP(i,j)=dW(i,j) + aE(i,j) + aN(i,j)  +fN +fE;
    aPv(j)=dW(i,j) + aE(i,j) + aN(i,j) + dS(i,j)  + fN +fE;

    %____________ CROSS DIFFUSION TERMS______________________

    %______________________U__________________________
            
    %Get velocity vertexes values
    

    %Split velocities
    u_wn=u_cor(i,j,1,1);
    u_en=u_cor(i,j,1,2);
    u_es=u_cor(i,j,1,3);
    u_ws=u_cor(i,j,1,4);

    %compute cross diffusion terms for each face 
    Sdu_w=dW_c(i,j)*(u_wn-u_ws);%face w
    Sdu_n=dN_c(i,j)*(u_en-u_wn);%face n
    Sdu_e=dE_c(i,j)*(u_en-u_es);%face e
    Sdu_s=dS_c(i,j)*(u_es-u_ws);%face s
    
    %Sum all contributions
    Sdu_T=Sdu_w+Sdu_n+Sdu_e+Sdu_s;

    %______________________V__________________________

    %Get velocity vertexes values
    

    %Split velocities
    v_wn=v_cor(i,j,1,1);
    v_en=v_cor(i,j,1,2);
    v_es=v_cor(i,j,1,3);
    v_ws=v_cor(i,j,1,4);

    %compute cross diffusion terms for each face 
    Sdv_w=dW_c(i,j)*(v_wn-v_ws);%face w
    Sdv_n=dN_c(i,j)*(v_en-v_wn);%face n
    Sdv_e=dE_c(i,j)*(v_en-v_es);%face e
    Sdv_s=dS_c(i,j)*(v_es-v_ws);%face s
    
    %Sum all contributions
    Sdv_T=Sdv_w+Sdv_n+Sdv_e+Sdv_s;

    %____________Get pressure Gradients__________________

    %get the gradient of pressure at the center of the cell
    
    
    %Grad in x
    grad_px=grad_p(i,j,1,1);
    %Grad in y
    grad_py=grad_p(i,j,1,2);

    %__________________TVD Corrections___________________________ 

    %Vector with mass fluxes 
    f_C=[fN,fE];

    %Vector with nodal values
    u_vec=[u(i-1,j),u(i,j+1),u(i,j)];
    v_vec=[v(i-1,j),v(i,j+1),v(i,j)];

    %vector with Gradients for U
    grad_u_mat=[grad_u(i-1,j,1,1),grad_u(i-1,j,1,2);...
        grad_u(i,j+1,1,1),grad_u(i,j+1,1,2);...
        grad_u(i,j,1,1),grad_u(i,j,1,2)];

    %vector with gradients for V
    grad_v_mat=[grad_v(i-1,j,1,1),grad_v(i-1,j,1,2);...
        grad_v(i,j+1,1,1),grad_v(i,j+1,1,2);...
        grad_v(i,j,1,1),grad_v(i,j,1,2)];

    %Vector with vectors from P to Neigborhood nodes
    r_mat=[VecCentNbNods(i,j,2,1),VecCentNbNods(i,j,2,2);...
        VecCentNbNods(i,j,3,1),VecCentNbNods(i,j,3,2)];

    [s_dc_total_x,s_dc_total_y] = tvdDefCorr(f_C,u_vec,v_vec,...
        grad_u_mat,grad_v_mat,r_mat);

    %Sources (Pressure dE(i,j)rivative)
    suX(i,j)=-cellVols(i,j)*grad_px + Sdu_T + u0*(dW(i,j) - fW) + s_dc_total_x;
    suY(i,j)=-cellVols(i,j)*grad_py + Sdv_T + s_dc_total_y;
        
end

