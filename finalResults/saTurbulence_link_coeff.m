function [Ant_W,Ant_N,Ant_E,Ant_S,Ant_P,Su_nt] = saTurbulence_link_coeff(...
    nu_tilde,nu,lower_omega,lgtFaces,cellVols,u_face,v_face,distMinWall,...
    grad_nu_tilde,dW,dN,dE,dS,dW_c,dN_c,dE_c,dS_c,nu_tilde_corners,...
    kappa,sigma_sa,cb1,cb2,cv1,cw1,cw2,cw3,...
    nx,ny,solidMask)
    % INPUTS
    %nu_tilde: nu_tilde from previous iteration
    %nu: Kinematic Viscosity
    %mu:Molecular viscosity
    %Rotation tensor norm 
    %lgtFaces: Area or Lenght (2D) of faces
    %u_face:velocity at faces w and e
    %v_face:velocity at faces n and s
    %grad_nu_tilde: gradient of nu tilde
    %,dW,dN,dE,dS :Direct diffusion coeffitients (already have mu)
    %dW_c,dN_c,dE_c,dS_c: Cross Diffusion coeffitiens (already have mu)
    %nu_tilde_corners: Nu tilde in the corners
    %kappa,sigma_sa,cb1,cb2,cv1,cw1,cw2,cw3: Spalart Allmaras model
    %constants
    %nx,ny: size of the matrices 
    %nxSolid,nxSymcond: Indexes for the solid part and free - stream zones

    %OUTPUTS
    %ant_W,ant_N,ant_E,ant_S,ant_P,suNt: Constants for linear systems of 
    %equations to solve for nu_tilde_k+1 

    %           %NU TILDE BOUNDARY CONDITIONS 
%                   NEUMMAN  DNU~/DY=0
%                                                      N  
%                                                      E   
%                                                      U  
% NU~=0.1NU                                            M  NU~/DX=0       
%                                                      M  
%                                                      A 
% ___ NU~=0.1NU ___|___   0    ___|___ NU~=0.1NU  ____ N

    %COEFFITIENTS FOR SPALART - ALLMARAS TURBULENCE MODEL EQUATION
    %Intialize variables
    Ant_W=zeros(ny,nx);
    Ant_N=zeros(ny,nx);
    Ant_E=zeros(ny,nx);
    Ant_S=zeros(ny,nx);
    Ant_P=zeros(ny,nx);
    Su_nt=zeros(ny,nx);

    %Elevate constant
    cw3_6=cw3^6;
    cv1_3=cv1^3;
    %_____________________________________________________________________
    %Compute total viscosity for difussion equation.

    nu_total=(nu+ nu_tilde)./sigma_sa;
    
    dW=dW.*(nu_total);
    dE=dE.*(nu_total);
    dN=dN.*(nu_total);
    dS=dS.*(nu_total);
    dW_c=dW_c.*(nu_total);
    dN_c=dN_c.*(nu_total);
    dE_c=dE_c.*(nu_total);
    dS_c=dS_c.*(nu_total);


    %1.- SPALART-ALLMARAS TURBULENCE NU TILDE TRANSPORT EQUATION
    %  LINK COEFFITIENTS
    
    %1.1- ################# INTERIOR CELLS ##############################
    for i=2:ny-1
        for j=2:nx-1

            %DIFFUSION AND CONVECTION
            %get face areas
            lgt_fw=lgtFaces(i,j,1);
            lgt_fn=lgtFaces(i,j,2); 
            lgt_fe=lgtFaces(i,j,3); 
            lgt_fs=lgtFaces(i,j,4);

            nu_tilde_k=nu_tilde(i,j); %Nu tilde from previous iteration

            %Volumetric fluxes at faces 
            fW=-lgt_fw*u_face(i,j);
            fN=-lgt_fn*v_face(i,j); 
            fE= lgt_fe*u_face(i,j+1);  
            fS= lgt_fs*v_face(i+1,j);

  
            %CROSS DIFFUSION EXPLICIT COMPUTATION
            %Get nu_tilde vertexes values

            %Split nu_tildes
            nu_tilde_wn=nu_tilde_corners(i,j,1,1);
            nu_tilde_en=nu_tilde_corners(i,j,1,2);
            nu_tilde_es=nu_tilde_corners(i,j,1,3);
            nu_tilde_ws=nu_tilde_corners(i,j,1,4);

            %compute cross diffusion terms for each face 
            Sdnu_t_w=dW_c(i,j)*(nu_tilde_wn-nu_tilde_ws);%face w
            Sdnu_t_n=dN_c(i,j)*(nu_tilde_en-nu_tilde_wn);%face n
            Sdnu_t_e=dE_c(i,j)*(nu_tilde_en-nu_tilde_es);%face e
            Sdnu_t_s=dS_c(i,j)*(nu_tilde_es-nu_tilde_ws);%face s

            %Sum all contributions
            Sdnu_t_T= Sdnu_t_w+ Sdnu_t_n + Sdnu_t_e + Sdnu_t_s;

            %Nonlinear Diffusion 
            dnu_t_dx=grad_nu_tilde(i,j,1,1);
            dnu_t_dy=grad_nu_tilde(i,j,1,2);
          
            nonLinDiff=(cb2/sigma_sa)*(dnu_t_dx^2 + dnu_t_dy^2);

            %LINEARIZATION OF SOURCE TERMS FOR IMPLICIT APROXIMATION
            %  Implicit Taylor 1st order expansion
            %Q(nu_tilde_k+1)=~ Q(nu_tilde_k) + Q'(nu_tilde_k)*(nu_tilde_k+1
            % - nu_tilde_k)
            %central coeffitient
            %AP=-Q'(nu_tilde_k)*cellVol
            %Source
            %AS=(Q(nu_tilde_k) - Q'(nu_tilde_k)*nu_tilde_k))*cellVol

            %Q(nu_tilde)=Tp(nu_tilde) + Td(nu_tilde)

            d_wall=distMinWall(i,j);

            %__________ Turbulence Production Tp(nu_tilde) ___________

            %Viscous ratio*
            x_nu=nu_tilde_k/nu;
            x_nu_prime=1/nu;
            
            %%Viscous damping function 1 fv1*
            fv1=(x_nu^3)/(x_nu^3 + cv1_3);
            fv1_prime=(3*x_nu_prime*(x_nu^2)*cv1_3)/((x_nu^3 + cv1_3)^2);

            %Viscous damping function 2 fv2*
            fv2=1 - x_nu/(1 + x_nu*fv1);
            fv2_prime=(((x_nu^2)*fv1_prime) - x_nu_prime)/((1 + x_nu*fv1)^2);
             
            %abs  vorticity (gamma tilde)*
            s_vort=lower_omega(i,j);
            
            %modified vorticity*           
            s_vort_tilde = s_vort + (nu_tilde_k*fv2)/((kappa*d_wall)^2);
            s_vort_tilde_prime= (1/((kappa*d_wall)^2))*(fv2 +...
                fv2_prime*nu_tilde_k);
            
            %tp Turbulence production*
            tp=cb1*s_vort_tilde*nu_tilde_k;
            tp_prime= cb1*(s_vort_tilde_prime*nu_tilde_k + s_vort_tilde);
            
            %__________ Turbulence Destruction Td(nu_tilde) __________
            
            %function r*
            r=min(10,nu_tilde_k/(s_vort_tilde * (kappa*d_wall)^2 ));
            r_prime=(1/((s_vort_tilde*kappa*d_wall)^2))*(s_vort_tilde -...
                nu_tilde_k*s_vort_tilde_prime);
            
            %function g*
            g=r + cw2*(r^6 -r);
            g_prime=r_prime*(1-cw2 + 6*cw2*r^5);
            
            %wall damping function fw
            fw=g*(((1 + cw3_6)/(g^6 + cw3_6))^(1/6));
            fw_prime=g_prime*cw3_6*(((1 + cw3_6)/((g^6 + cw3_6)^7))^(1/6));

            %Td
            td=-cw1*fw*(nu_tilde_k/d_wall)^2;
            td_prime=-cw1*(fw_prime*(nu_tilde_k/d_wall)^2 + ...
                2*fw*(nu_tilde_k/(d_wall^2)));
            
            %ALL SOURCES Q(nu_tilde_k) & Q'(nu_tilde_k)
            q_s=tp + td;
            q_s_prime=tp_prime + td_prime;

            %Coeffitiens using upwind squeme

            Ant_W(i,j)=dW(i,j) + max(0,-fW);
            Ant_N(i,j)=dN(i,j) + max(0,-fE);
            Ant_E(i,j)=dE(i,j) + max(0,-fN);
            Ant_S(i,j)=dS(i,j) + max(0,-fS);

            Ant_P(i,j)=Ant_W(i,j) + Ant_N(i,j) + Ant_E(i,j) +Ant_S(i,j) ...
                +fW +fN +fE +fS - cellVols(i,j)*q_s_prime;

            Su_nt(i,j)=Sdnu_t_T + cellVols(i,j)*(nonLinDiff + q_s - ...
                q_s_prime*nu_tilde_k);

        end
    end

    %#########################   EDGES ##################################
    %West Edge_________________________________________________________
    %FIXED NU TILDE NU~=0.1*NU
    j=1;
    for i=2:ny-1

        %DIFFUSION AND CONVECTION
        %get face areas
        lgt_fw=lgtFaces(i,j,1);
        lgt_fn=lgtFaces(i,j,2); 
        lgt_fe=lgtFaces(i,j,3); 
        lgt_fs=lgtFaces(i,j,4);

        nu_tilde_k=nu_tilde(i,j); %Nu tilde from previous iteration

       
        %Volumetric fluxes at faces 
        fW=-lgt_fw*u_face(i,j);
        fN=-lgt_fn*v_face(i,j); 
        fE= lgt_fe*u_face(i,j+1);       
        fS= lgt_fs*v_face(i+1,j);

        %CROSS DIFFUSION EXPLICIT COMPUTATION
        %Get nu_tilde vertexes values
        
        %Split nu_tildes
        nu_tilde_wn=nu_tilde_corners(i,j,1,1);
        nu_tilde_en=nu_tilde_corners(i,j,1,2);
        nu_tilde_es=nu_tilde_corners(i,j,1,3);
        nu_tilde_ws=nu_tilde_corners(i,j,1,4);

        %compute cross diffusion terms for each face 
        Sdnu_t_w=dW_c(i,j)*(nu_tilde_wn-nu_tilde_ws);%face w
        Sdnu_t_n=dN_c(i,j)*(nu_tilde_en-nu_tilde_wn);%face n
        Sdnu_t_e=dE_c(i,j)*(nu_tilde_en-nu_tilde_es);%face e
        Sdnu_t_s=dS_c(i,j)*(nu_tilde_es-nu_tilde_ws);%face s

        %Sum all contributions
        Sdnu_t_T= Sdnu_t_w+ Sdnu_t_n + Sdnu_t_e + Sdnu_t_s;

        %Nonlinear Diffusion 
        dnu_t_dx=grad_nu_tilde(i,j,1,1);
        dnu_t_dy=grad_nu_tilde(i,j,1,2);
      
        nonLinDiff=(cb2/sigma_sa)*(dnu_t_dx^2 + dnu_t_dy^2);

        %LINEARIZATION OF SOURCE TERMS FOR IMPLICIT APROXIMATION
        %  Implicit Taylor 1st order expansion
        %Q(nu_tilde_k+1)=~ Q(nu_tilde_k) + Q'(nu_tilde_k)*(nu_tilde_k+1
        % - nu_tilde_k)
        %central coeffitient
        %AP=-Q'(nu_tilde_k)*cellVol
        %Source
        %AS=(Q(nu_tilde_k) - Q'(nu_tilde_k)*nu_tilde_k))*cellVol

        %Q(nu_tilde)=Tp(nu_tilde) + Td(nu_tilde)

        d_wall=distMinWall(i,j);

        %__________ Turbulence Production Tp(nu_tilde) ___________

        %Viscous ratio*
        x_nu=nu_tilde_k/nu;
        x_nu_prime=1/nu;
        
        %%Viscous damping function 1 fv1*
        fv1=(x_nu^3)/(x_nu^3 + cv1_3);
        fv1_prime=(3*x_nu_prime*(x_nu^2)*cv1_3)/((x_nu^3 + cv1_3)^2);

        %Viscous damping function 2 fv2*
        fv2=1 - x_nu/(1 + x_nu*fv1);
        fv2_prime=(((x_nu^2)*fv1_prime) - x_nu_prime)/((1 + x_nu*fv1)^2);
         
        %abs  vorticity (gamma tilde)*
        s_vort=lower_omega(i,j);
        
        %modified vorticity*           
        s_vort_tilde = s_vort + (nu_tilde_k*fv2)/((kappa*d_wall)^2);
        s_vort_tilde_prime= (1/((kappa*d_wall)^2))*(fv2 +...
            fv2_prime*nu_tilde_k);
        
        %tp Turbulence production*
        tp=cb1*s_vort_tilde*nu_tilde_k;
        tp_prime= cb1*(s_vort_tilde_prime*nu_tilde_k + s_vort_tilde);
        
        %__________ Turbulence Destruction Td(nu_tilde) __________
        
        %function r*
        r=min(10,nu_tilde_k/(s_vort_tilde * (kappa*d_wall)^2 ));
        r_prime=(1/((s_vort_tilde*kappa*d_wall)^2))*(s_vort_tilde -...
            nu_tilde_k*s_vort_tilde_prime);
        
        %function g*
        g=r + cw2*(r^6 -r);
        g_prime=r_prime*(1-cw2 + 6*cw2*r^5);
        
        %wall damping function fw
        fw=g*(((1 + cw3_6)/(g^6 + cw3_6))^(1/6));
        fw_prime=g_prime*cw3_6*(((1 + cw3_6)/((g^6 + cw3_6)^7))^(1/6));

        %Td
        td=-cw1*fw*(nu_tilde_k/d_wall)^2;
        td_prime=-cw1*(fw_prime*(nu_tilde_k/d_wall)^2 + ...
            2*fw*(nu_tilde_k/(d_wall^2)));
        
        %ALL SOURCES Q(nu_tilde_k) & Q'(nu_tilde_k)
        q_s=tp + td;
        q_s_prime=tp_prime + td_prime;

        %Coeffitiens using upwind squeme

        Ant_W(i,j)=0;
        Ant_N(i,j)=dN(i,j) + max(0,-fE);
        Ant_E(i,j)=dE(i,j) + max(0,-fN);
        Ant_S(i,j)=dS(i,j) + max(0,-fS);

        Ant_P(i,j)=dW(i,j) + Ant_N(i,j) + Ant_E(i,j) +Ant_S(i,j) ...
             +fN +fE +fS - cellVols(i,j)*q_s_prime;

        Su_nt(i,j)=Sdnu_t_T + cellVols(i,j)*(nonLinDiff + q_s - ...
            q_s_prime*nu_tilde_k) + 0.1*nu*(dW(i,j)- fW);

    end

    %North Edge_________________________________________________________
    %NEUMMAN
    i=1;
    for j=2:nx-1
        
        %DIFFUSION AND CONVECTION
        %get face areas
        lgt_fw=lgtFaces(i,j,1);
        lgt_fn=lgtFaces(i,j,2); 
        lgt_fe=lgtFaces(i,j,3); 
        lgt_fs=lgtFaces(i,j,4);

        nu_tilde_k=nu_tilde(i,j); %Nu tilde from previous iteration

        %Volumetric fluxes at faces 
        fW=-lgt_fw*u_face(i,j);
        fN=-lgt_fn*v_face(i,j); 
        fE= lgt_fe*u_face(i,j+1); 
        fS= lgt_fs*v_face(i+1,j);

        %CROSS DIFFUSION EXPLICIT COMPUTATION
        %Get nu_tilde vertexes values
        
        %Split nu_tildes
        nu_tilde_wn=nu_tilde_corners(i,j,1,1);
        nu_tilde_en=nu_tilde_corners(i,j,1,2);
        nu_tilde_es=nu_tilde_corners(i,j,1,3);
        nu_tilde_ws=nu_tilde_corners(i,j,1,4);

        %compute cross diffusion terms for each face 
        Sdnu_t_w=dW_c(i,j)*(nu_tilde_wn-nu_tilde_ws);%face w
        Sdnu_t_n=dN_c(i,j)*(nu_tilde_en-nu_tilde_wn);%face n
        Sdnu_t_e=dE_c(i,j)*(nu_tilde_en-nu_tilde_es);%face e
        Sdnu_t_s=dS_c(i,j)*(nu_tilde_es-nu_tilde_ws);%face s

        %Sum all contributions
        Sdnu_t_T= Sdnu_t_w+ Sdnu_t_n + Sdnu_t_e + Sdnu_t_s;

        %Nonlinear Diffusion 
        dnu_t_dx=grad_nu_tilde(i,j,1,1);
        dnu_t_dy=grad_nu_tilde(i,j,1,2);
      
        nonLinDiff=(cb2/sigma_sa)*(dnu_t_dx^2 + dnu_t_dy^2);

        %LINEARIZATION OF SOURCE TERMS FOR IMPLICIT APROXIMATION
        %  Implicit Taylor 1st order expansion
        %Q(nu_tilde_k+1)=~ Q(nu_tilde_k) + Q'(nu_tilde_k)*(nu_tilde_k+1
        % - nu_tilde_k)
        %central coeffitient
        %AP=-Q'(nu_tilde_k)*cellVol
        %Source
        %AS=(Q(nu_tilde_k) - Q'(nu_tilde_k)*nu_tilde_k))*cellVol

        %Q(nu_tilde)=Tp(nu_tilde) + Td(nu_tilde)

        d_wall=distMinWall(i,j);

        %__________ Turbulence Production Tp(nu_tilde) ___________

        %Viscous ratio*
        x_nu=nu_tilde_k/nu;
        x_nu_prime=1/nu;
        
        %%Viscous damping function 1 fv1*
        fv1=(x_nu^3)/(x_nu^3 + cv1_3);
        fv1_prime=(3*x_nu_prime*(x_nu^2)*cv1_3)/((x_nu^3 + cv1_3)^2);

        %Viscous damping function 2 fv2*
        fv2=1 - x_nu/(1 + x_nu*fv1);
        fv2_prime=(((x_nu^2)*fv1_prime) - x_nu_prime)/((1 + x_nu*fv1)^2);
         
        %abs  vorticity (gamma tilde)*
        s_vort=lower_omega(i,j);
        
        %modified vorticity*           
        s_vort_tilde = s_vort + (nu_tilde_k*fv2)/((kappa*d_wall)^2);
        s_vort_tilde_prime= (1/((kappa*d_wall)^2))*(fv2 +...
            fv2_prime*nu_tilde_k);
        
        %tp Turbulence production*
        tp=cb1*s_vort_tilde*nu_tilde_k;
        tp_prime= cb1*(s_vort_tilde_prime*nu_tilde_k + s_vort_tilde);
        
        %__________ Turbulence Destruction Td(nu_tilde) __________
        
        %function r*
        r=min(10,nu_tilde_k/(s_vort_tilde * (kappa*d_wall)^2 ));
        r_prime=(1/((s_vort_tilde*kappa*d_wall)^2))*(s_vort_tilde -...
            nu_tilde_k*s_vort_tilde_prime);
        
        %function g*
        g=r + cw2*(r^6 -r);
        g_prime=r_prime*(1-cw2 + 6*cw2*r^5);
        
        %wall damping function fw
        fw=g*(((1 + cw3_6)/(g^6 + cw3_6))^(1/6));
        fw_prime=g_prime*cw3_6*(((1 + cw3_6)/((g^6 + cw3_6)^7))^(1/6));

        %Td
        td=-cw1*fw*(nu_tilde_k/d_wall)^2;
        td_prime=-cw1*(fw_prime*(nu_tilde_k/d_wall)^2 + ...
            2*fw*(nu_tilde_k/(d_wall^2)));
        
        %ALL SOURCES Q(nu_tilde_k) & Q'(nu_tilde_k)
        q_s=tp + td;
        q_s_prime=tp_prime + td_prime;

        %Coeffitiens using upwind squeme

        Ant_W(i,j)=dW(i,j) + max(0,-fW);
        Ant_N(i,j)=0;
        Ant_E(i,j)=dE(i,j) + max(0,-fN);
        Ant_S(i,j)=dS(i,j) + max(0,-fS);

        Ant_P(i,j)=Ant_W(i,j) + Ant_E(i,j) +Ant_S(i,j) ...
            +fW +fN +fE +fS - cellVols(i,j)*q_s_prime;

        Su_nt(i,j)=Sdnu_t_T + cellVols(i,j)*(nonLinDiff + q_s - ...
            q_s_prime*nu_tilde_k);

    end
    %East Edge_________________________________________________________
    %NEUMMAN
    j=nx;
    for i=2:ny-1
        
        %DIFFUSION AND CONVECTION
        %get face areas
        lgt_fw=lgtFaces(i,j,1);
        lgt_fn=lgtFaces(i,j,2); 
        lgt_fe=lgtFaces(i,j,3); 
        lgt_fs=lgtFaces(i,j,4);
    
        nu_tilde_k=nu_tilde(i,j); %Nu tilde from previous iteration

        %Volumetric fluxes at faces 
        fW=-lgt_fw*u_face(i,j);
        fN=-lgt_fn*v_face(i,j);
        fE= lgt_fe*u_face(i,j+1);          
        fS= lgt_fs*v_face(i+1,j);
        
        %CROSS DIFFUSION EXPLICIT COMPUTATION
        %Get nu_tilde vertexes values        
    
        %Split nu_tildes
        nu_tilde_wn=nu_tilde_corners(i,j,1,1);
        nu_tilde_en=nu_tilde_corners(i,j,1,2);
        nu_tilde_es=nu_tilde_corners(i,j,1,3);
        nu_tilde_ws=nu_tilde_corners(i,j,1,4);
    
        %compute cross diffusion terms for each face 
        Sdnu_t_w=dW_c(i,j)*(nu_tilde_wn-nu_tilde_ws);%face w
        Sdnu_t_n=dN_c(i,j)*(nu_tilde_en-nu_tilde_wn);%face n
        Sdnu_t_e=dE_c(i,j)*(nu_tilde_en-nu_tilde_es);%face e
        Sdnu_t_s=dS_c(i,j)*(nu_tilde_es-nu_tilde_ws);%face s
    
        %Sum all contributions
        Sdnu_t_T= Sdnu_t_w+ Sdnu_t_n + Sdnu_t_e + Sdnu_t_s;
    
        %Nonlinear Diffusion 
        dnu_t_dx=grad_nu_tilde(i,j,1,1);
        dnu_t_dy=grad_nu_tilde(i,j,1,2);
      
        nonLinDiff=(cb2/sigma_sa)*(dnu_t_dx^2 + dnu_t_dy^2);
    
        %LINEARIZATION OF SOURCE TERMS FOR IMPLICIT APROXIMATION
        %  Implicit Taylor 1st order expansion
        %Q(nu_tilde_k+1)=~ Q(nu_tilde_k) + Q'(nu_tilde_k)*(nu_tilde_k+1
        % - nu_tilde_k)
        %central coeffitient
        %AP=-Q'(nu_tilde_k)*cellVol
        %Source
        %AS=(Q(nu_tilde_k) - Q'(nu_tilde_k)*nu_tilde_k))*cellVol

        %Q(nu_tilde)=Tp(nu_tilde) + Td(nu_tilde)

        d_wall=distMinWall(i,j);

        %__________ Turbulence Production Tp(nu_tilde) ___________

        %Viscous ratio*
        x_nu=nu_tilde_k/nu;
        x_nu_prime=1/nu;
        
        %%Viscous damping function 1 fv1*
        fv1=(x_nu^3)/(x_nu^3 + cv1_3);
        fv1_prime=(3*x_nu_prime*(x_nu^2)*cv1_3)/((x_nu^3 + cv1_3)^2);

        %Viscous damping function 2 fv2*
        fv2=1 - x_nu/(1 + x_nu*fv1);
        fv2_prime=(((x_nu^2)*fv1_prime) - x_nu_prime)/((1 + x_nu*fv1)^2);
         
        %abs  vorticity (gamma tilde)*
        s_vort=lower_omega(i,j);
        
        %modified vorticity*           
        s_vort_tilde = s_vort + (nu_tilde_k*fv2)/((kappa*d_wall)^2);
        s_vort_tilde_prime= (1/((kappa*d_wall)^2))*(fv2 +...
            fv2_prime*nu_tilde_k);
        
        %tp Turbulence production*
        tp=cb1*s_vort_tilde*nu_tilde_k;
        tp_prime= cb1*(s_vort_tilde_prime*nu_tilde_k + s_vort_tilde);
        
        %__________ Turbulence Destruction Td(nu_tilde) __________
        
        %function r*
        r=min(10,nu_tilde_k/(s_vort_tilde * (kappa*d_wall)^2 ));
        r_prime=(1/((s_vort_tilde*kappa*d_wall)^2))*(s_vort_tilde -...
            nu_tilde_k*s_vort_tilde_prime);
        
        %function g*
        g=r + cw2*(r^6 -r);
        g_prime=r_prime*(1-cw2 + 6*cw2*r^5);
        
        %wall damping function fw
        fw=g*(((1 + cw3_6)/(g^6 + cw3_6))^(1/6));
        fw_prime=g_prime*cw3_6*(((1 + cw3_6)/((g^6 + cw3_6)^7))^(1/6));

        %Td
        td=-cw1*fw*(nu_tilde_k/d_wall)^2;
        td_prime=-cw1*(fw_prime*(nu_tilde_k/d_wall)^2 + ...
            2*fw*(nu_tilde_k/(d_wall^2)));
        
        %ALL SOURCES Q(nu_tilde_k) & Q'(nu_tilde_k)
        q_s=tp + td;
        q_s_prime=tp_prime + td_prime;
    
        %Coeffitiens using upwind squeme
    
        Ant_W(i,j)=dW(i,j) + max(0,-fW);
        Ant_N(i,j)=dN(i,j) + max(0,-fE);
        Ant_E(i,j)=0;
        Ant_S(i,j)=dS(i,j) + max(0,-fS);
    
        Ant_P(i,j)=Ant_W(i,j) + Ant_N(i,j)  +Ant_S(i,j) ...
            +fW +fN +fE +fS - cellVols(i,j)*q_s_prime;
    
        Su_nt(i,j)=Sdnu_t_T + cellVols(i,j)*(nonLinDiff + q_s - ...
            q_s_prime*nu_tilde_k);

    end
    %South Edge_________________________________________________________
    %FIXED NU~ AT FREE STREAM 0 AT SOLID SURFACE
    i=ny;
    for j=2:nx-1

        %get face areas
        lgt_fw=lgtFaces(i,j,1);
        lgt_fn=lgtFaces(i,j,2); 
        lgt_fe=lgtFaces(i,j,3); 
        lgt_fs=lgtFaces(i,j,4);

        nu_tilde_k=nu_tilde(i,j); %Nu tilde from previous iteration

        %Volumetric fluxes at faces 
        fW=-lgt_fw*u_face(i,j);
        fN=-lgt_fn*v_face(i,j);
        fE= lgt_fe*u_face(i,j+1);         
        fS=-lgt_fs*v_face(i+1,j);

        %CROSS DIFFUSION EXPLICIT COMPUTATION
        %Get nu_tilde vertexes values
        
        %Split nu_tildes
        nu_tilde_wn=nu_tilde_corners(i,j,1,1);
        nu_tilde_en=nu_tilde_corners(i,j,1,2);
        nu_tilde_es=nu_tilde_corners(i,j,1,3);
        nu_tilde_ws=nu_tilde_corners(i,j,1,4);

        %compute cross diffusion terms for each face 
        Sdnu_t_w=dW_c(i,j)*(nu_tilde_wn-nu_tilde_ws);%face w
        Sdnu_t_n=dN_c(i,j)*(nu_tilde_en-nu_tilde_wn);%face n
        Sdnu_t_e=dE_c(i,j)*(nu_tilde_en-nu_tilde_es);%face e
        Sdnu_t_s=dS_c(i,j)*(nu_tilde_es-nu_tilde_ws);%face s

        %Sum all contributions
        Sdnu_t_T= Sdnu_t_w+ Sdnu_t_n + Sdnu_t_e + Sdnu_t_s;

        %Nonlinear Diffusion 
        dnu_t_dx=grad_nu_tilde(i,j,1,1);
        dnu_t_dy=grad_nu_tilde(i,j,1,2);
      
        nonLinDiff=(cb2/sigma_sa)*(dnu_t_dx^2 + dnu_t_dy^2);

        %LINEARIZATION OF SOURCE TERMS FOR IMPLICIT APROXIMATION
        %  Implicit Taylor 1st order expansion
        %Q(nu_tilde_k+1)=~ Q(nu_tilde_k) + Q'(nu_tilde_k)*(nu_tilde_k+1
        % - nu_tilde_k)
        %central coeffitient
        %AP=-Q'(nu_tilde_k)*cellVol
        %Source
        %AS=(Q(nu_tilde_k) - Q'(nu_tilde_k)*nu_tilde_k))*cellVol

        %Q(nu_tilde)=Tp(nu_tilde) + Td(nu_tilde)

        d_wall=distMinWall(i,j);

        %__________ Turbulence Production Tp(nu_tilde) ___________

        %Viscous ratio*
        x_nu=nu_tilde_k/nu;
        x_nu_prime=1/nu;
        
        %%Viscous damping function 1 fv1*
        fv1=(x_nu^3)/(x_nu^3 + cv1_3);
        fv1_prime=(3*x_nu_prime*(x_nu^2)*cv1_3)/((x_nu^3 + cv1_3)^2);

        %Viscous damping function 2 fv2*
        fv2=1 - x_nu/(1 + x_nu*fv1);
        fv2_prime=(((x_nu^2)*fv1_prime) - x_nu_prime)/((1 + x_nu*fv1)^2);
         
        %abs  vorticity (gamma tilde)*
        s_vort=lower_omega(i,j);
        
        %modified vorticity*           
        s_vort_tilde = s_vort + (nu_tilde_k*fv2)/((kappa*d_wall)^2);
        s_vort_tilde_prime= (1/((kappa*d_wall)^2))*(fv2 +...
            fv2_prime*nu_tilde_k);
        
        %tp Turbulence production*
        tp=cb1*s_vort_tilde*nu_tilde_k;
        tp_prime= cb1*(s_vort_tilde_prime*nu_tilde_k + s_vort_tilde);
        
        %__________ Turbulence Destruction Td(nu_tilde) __________
        
        %function r*
        r=min(10,nu_tilde_k/(s_vort_tilde * (kappa*d_wall)^2 ));
        r_prime=(1/((s_vort_tilde*kappa*d_wall)^2))*(s_vort_tilde -...
            nu_tilde_k*s_vort_tilde_prime);
        
        %function g*
        g=r + cw2*(r^6 -r);
        g_prime=r_prime*(1-cw2 + 6*cw2*r^5);
        
        %wall damping function fw
        fw=g*(((1 + cw3_6)/(g^6 + cw3_6))^(1/6));
        fw_prime=g_prime*cw3_6*(((1 + cw3_6)/((g^6 + cw3_6)^7))^(1/6));

        %Td
        td=-cw1*fw*(nu_tilde_k/d_wall)^2;
        td_prime=-cw1*(fw_prime*(nu_tilde_k/d_wall)^2 + ...
            2*fw*(nu_tilde_k/(d_wall^2)));
        
        %ALL SOURCES Q(nu_tilde_k) & Q'(nu_tilde_k)
        q_s=tp + td;
        q_s_prime=tp_prime + td_prime;

        %Coeffitiens using upwind squeme

        Ant_W(i,j)=dW(i,j) + max(0,-fW);
        Ant_N(i,j)=dN(i,j) + max(0,-fE);
        Ant_E(i,j)=dE(i,j) + max(0,-fN);
        Ant_S(i,j)=0;

        Su_nt(i,j)=Sdnu_t_T + cellVols(i,j)*(nonLinDiff + q_s - ...
            q_s_prime*nu_tilde_k);

        if solidMask(j)% solid surf.

            Ant_P(i,j)=Ant_W(i,j) + Ant_N(i,j) + Ant_E(i,j) +dS(i,j) ...
                +fW +fN +fE - cellVols(i,j)*q_s_prime;

            Su_nt(i,j)=Sdnu_t_T + cellVols(i,j)*(nonLinDiff + q_s - ...
                q_s_prime*nu_tilde_k);

        else %free stream

            Ant_P(i,j)=Ant_W(i,j) + Ant_N(i,j) + Ant_E(i,j) + dS(i,j)...
                +fW +fN +fE - cellVols(i,j)*q_s_prime;

            Su_nt(i,j)=Sdnu_t_T + cellVols(i,j)*(nonLinDiff + q_s - ...
                q_s_prime*nu_tilde_k) + 0.1*nu*(dS(i,j)-fS);

        end

    end
    %########################## Corners ##################################
    %west north_________________________________________________________
    %FIXED AND NEUMMAN
    i=1;
    j=1;
    %DIFFUSION AND CONVECTION
    %get face areas
    lgt_fw=lgtFaces(i,j,1);
    lgt_fn=lgtFaces(i,j,2); 
    lgt_fe=lgtFaces(i,j,3); 
    lgt_fs=lgtFaces(i,j,4);

    nu_tilde_k=nu_tilde(i,j); %Nu tilde from previous iteration

    %Volumetric fluxes at faces 
    fW=-lgt_fw*u_face(i,j);
    fN=-lgt_fn*v_face(i,j); 
    fE= lgt_fe*u_face(i,j+1); 
    fS= lgt_fs*v_face(i+1,j);

    %CROSS DIFFUSION EXPLICIT COMPUTATION
    %Get nu_tilde vertexes values
    
    %Split nu_tildes
    nu_tilde_wn=nu_tilde_corners(i,j,1,1);
    nu_tilde_en=nu_tilde_corners(i,j,1,2);
    nu_tilde_es=nu_tilde_corners(i,j,1,3);
    nu_tilde_ws=nu_tilde_corners(i,j,1,4);

    %compute cross diffusion terms for each face 
    Sdnu_t_w=dW_c(i,j)*(nu_tilde_wn-nu_tilde_ws);%face w
    Sdnu_t_n=dN_c(i,j)*(nu_tilde_en-nu_tilde_wn);%face n
    Sdnu_t_e=dE_c(i,j)*(nu_tilde_en-nu_tilde_es);%face e
    Sdnu_t_s=dS_c(i,j)*(nu_tilde_es-nu_tilde_ws);%face s

    %Sum all contributions
    Sdnu_t_T= Sdnu_t_w+ Sdnu_t_n + Sdnu_t_e + Sdnu_t_s;

    %Nonlinear Diffusion 
    dnu_t_dx=grad_nu_tilde(i,j,1,1);
    dnu_t_dy=grad_nu_tilde(i,j,1,2);
  
    nonLinDiff=(cb2/sigma_sa)*(dnu_t_dx^2 + dnu_t_dy^2);

    %LINEARIZATION OF SOURCE TERMS FOR IMPLICIT APROXIMATION
    %  Implicit Taylor 1st order expansion
    %Q(nu_tilde_k+1)=~ Q(nu_tilde_k) + Q'(nu_tilde_k)*(nu_tilde_k+1
    % - nu_tilde_k)
    %central coeffitient
    %AP=-Q'(nu_tilde_k)*cellVol
    %Source
    %AS=(Q(nu_tilde_k) - Q'(nu_tilde_k)*nu_tilde_k))*cellVol

    %Q(nu_tilde)=Tp(nu_tilde) + Td(nu_tilde)

    d_wall=distMinWall(i,j);

    %__________ Turbulence Production Tp(nu_tilde) ___________

    %Viscous ratio*
    x_nu=nu_tilde_k/nu;
    x_nu_prime=1/nu;
    
    %%Viscous damping function 1 fv1*
    fv1=(x_nu^3)/(x_nu^3 + cv1_3);
    fv1_prime=(3*x_nu_prime*(x_nu^2)*cv1_3)/((x_nu^3 + cv1_3)^2);

    %Viscous damping function 2 fv2*
    fv2=1 - x_nu/(1 + x_nu*fv1);
    fv2_prime=(((x_nu^2)*fv1_prime) - x_nu_prime)/((1 + x_nu*fv1)^2);
     
    %abs  vorticity (gamma tilde)*
    s_vort=lower_omega(i,j);
    
    %modified vorticity*           
    s_vort_tilde = s_vort + (nu_tilde_k*fv2)/((kappa*d_wall)^2);
    s_vort_tilde_prime= (1/((kappa*d_wall)^2))*(fv2 +...
        fv2_prime*nu_tilde_k);
    
    %tp Turbulence production*
    tp=cb1*s_vort_tilde*nu_tilde_k;
    tp_prime= cb1*(s_vort_tilde_prime*nu_tilde_k + s_vort_tilde);
    
    %__________ Turbulence Destruction Td(nu_tilde) __________
    
    %function r*
    r=min(10,nu_tilde_k/(s_vort_tilde * (kappa*d_wall)^2 ));
    r_prime=(1/((s_vort_tilde*kappa*d_wall)^2))*(s_vort_tilde -...
        nu_tilde_k*s_vort_tilde_prime);
    
    %function g*
    g=r + cw2*(r^6 -r);
    g_prime=r_prime*(1-cw2 + 6*cw2*r^5);
    
    %wall damping function fw
    fw=g*(((1 + cw3_6)/(g^6 + cw3_6))^(1/6));
    fw_prime=g_prime*cw3_6*(((1 + cw3_6)/((g^6 + cw3_6)^7))^(1/6));

    %Td
    td=-cw1*fw*(nu_tilde_k/d_wall)^2;
    td_prime=-cw1*(fw_prime*(nu_tilde_k/d_wall)^2 + ...
        2*fw*(nu_tilde_k/(d_wall^2)));
    
    %ALL SOURCES Q(nu_tilde_k) & Q'(nu_tilde_k)
    q_s=tp + td;
    q_s_prime=tp_prime + td_prime;

    %Coeffitiens using upwind squeme

    Ant_W(i,j)=0;
    Ant_N(i,j)=0;
    Ant_E(i,j)=dE(i,j) + max(0,-fN);
    Ant_S(i,j)=dS(i,j) + max(0,-fS);

    Ant_P(i,j)=dW(i,j) + Ant_E(i,j) +Ant_S(i,j) ...
        +fN +fE +fS - cellVols(i,j)*q_s_prime;

    Su_nt(i,j)=Sdnu_t_T + cellVols(i,j)*(nonLinDiff + q_s - ...
        q_s_prime*nu_tilde_k) + 0.1*nu*(dW(i,j)- fW) ;


    %east north_________________________________________________________
    %NEUMMAN X2
    i=1;
    j=nx;
    
    %DIFFUSION AND CONVECTION
    %get face areas
    lgt_fw=lgtFaces(i,j,1);
    lgt_fn=lgtFaces(i,j,2); 
    lgt_fe=lgtFaces(i,j,3); 
    lgt_fs=lgtFaces(i,j,4);

    nu_tilde_k=nu_tilde(i,j); %Nu tilde from previous iteration

    %Volumetric fluxes at faces 
    fW=-lgt_fw*u_face(i,j);
    fN=-lgt_fn*v_face(i,j); 
    fE= lgt_fe*u_face(i,j+1);    
    fS= lgt_fs*v_face(i+1,j);

    %CROSS DIFFUSION EXPLICIT COMPUTATION
    %Get nu_tilde vertexes values    

    %Split nu_tildes
    nu_tilde_wn=nu_tilde_corners(i,j,1,1);
    nu_tilde_en=nu_tilde_corners(i,j,1,2);
    nu_tilde_es=nu_tilde_corners(i,j,1,3);
    nu_tilde_ws=nu_tilde_corners(i,j,1,4);

    %compute cross diffusion terms for each face 
    Sdnu_t_w=dW_c(i,j)*(nu_tilde_wn-nu_tilde_ws);%face w
    Sdnu_t_n=dN_c(i,j)*(nu_tilde_en-nu_tilde_wn);%face n
    Sdnu_t_e=dE_c(i,j)*(nu_tilde_en-nu_tilde_es);%face e
    Sdnu_t_s=dS_c(i,j)*(nu_tilde_es-nu_tilde_ws);%face s

    %Sum all contributions
    Sdnu_t_T= Sdnu_t_w+ Sdnu_t_n + Sdnu_t_e + Sdnu_t_s;

    %Nonlinear Diffusion 
    dnu_t_dx=grad_nu_tilde(i,j,1,1);
    dnu_t_dy=grad_nu_tilde(i,j,1,2);
  
    nonLinDiff=(cb2/sigma_sa)*(dnu_t_dx^2 + dnu_t_dy^2);

    %LINEARIZATION OF SOURCE TERMS FOR IMPLICIT APROXIMATION
    %  Implicit Taylor 1st order expansion
    %Q(nu_tilde_k+1)=~ Q(nu_tilde_k) + Q'(nu_tilde_k)*(nu_tilde_k+1
    % - nu_tilde_k)
    %central coeffitient
    %AP=-Q'(nu_tilde_k)*cellVol
    %Source
    %AS=(Q(nu_tilde_k) - Q'(nu_tilde_k)*nu_tilde_k))*cellVol

    %Q(nu_tilde)=Tp(nu_tilde) + Td(nu_tilde)

    d_wall=distMinWall(i,j);

    %__________ Turbulence Production Tp(nu_tilde) ___________

    %Viscous ratio*
    x_nu=nu_tilde_k/nu;
    x_nu_prime=1/nu;
    
    %%Viscous damping function 1 fv1*
    fv1=(x_nu^3)/(x_nu^3 + cv1_3);
    fv1_prime=(3*x_nu_prime*(x_nu^2)*cv1_3)/((x_nu^3 + cv1_3)^2);

    %Viscous damping function 2 fv2*
    fv2=1 - x_nu/(1 + x_nu*fv1);
    fv2_prime=(((x_nu^2)*fv1_prime) - x_nu_prime)/((1 + x_nu*fv1)^2);
     
    %abs  vorticity (gamma tilde)*
    s_vort=lower_omega(i,j);
    
    %modified vorticity*           
    s_vort_tilde = s_vort + (nu_tilde_k*fv2)/((kappa*d_wall)^2);
    s_vort_tilde_prime= (1/((kappa*d_wall)^2))*(fv2 +...
        fv2_prime*nu_tilde_k);
    
    %tp Turbulence production*
    tp=cb1*s_vort_tilde*nu_tilde_k;
    tp_prime= cb1*(s_vort_tilde_prime*nu_tilde_k + s_vort_tilde);
    
    %__________ Turbulence Destruction Td(nu_tilde) __________
    
    %function r*
    r=min(10,nu_tilde_k/(s_vort_tilde * (kappa*d_wall)^2 ));
    r_prime=(1/((s_vort_tilde*kappa*d_wall)^2))*(s_vort_tilde -...
        nu_tilde_k*s_vort_tilde_prime);
    
    %function g*
    g=r + cw2*(r^6 -r);
    g_prime=r_prime*(1-cw2 + 6*cw2*r^5);
    
    %wall damping function fw
    fw=g*(((1 + cw3_6)/(g^6 + cw3_6))^(1/6));
    fw_prime=g_prime*cw3_6*(((1 + cw3_6)/((g^6 + cw3_6)^7))^(1/6));

    %Td
    td=-cw1*fw*(nu_tilde_k/d_wall)^2;
    td_prime=-cw1*(fw_prime*(nu_tilde_k/d_wall)^2 + ...
        2*fw*(nu_tilde_k/(d_wall^2)));
    
    %ALL SOURCES Q(nu_tilde_k) & Q'(nu_tilde_k)
    q_s=tp + td;
    q_s_prime=tp_prime + td_prime;

    %Coeffitiens using upwind squeme

    Ant_W(i,j)=dW(i,j) + max(0,-fW);
    Ant_N(i,j)=0;
    Ant_E(i,j)=0;
    Ant_S(i,j)=dS(i,j) + max(0,-fS);

    Ant_P(i,j)=Ant_W(i,j) +Ant_S(i,j)+fW +fN +fE +fS - ...
        cellVols(i,j)*q_s_prime;

    Su_nt(i,j)=Sdnu_t_T + cellVols(i,j)*(nonLinDiff + q_s - ...
        q_s_prime*nu_tilde_k);

    %east south_________________________________________________________
    %NEUMMAN AND FIXED
    i=ny;
    j=nx;

    %DIFFUSION AND CONVECTION
    %get face areas
    lgt_fw=lgtFaces(i,j,1);
    lgt_fn=lgtFaces(i,j,2); 
    lgt_fe=lgtFaces(i,j,3); 
    lgt_fs=lgtFaces(i,j,4);

    nu_tilde_k=nu_tilde(i,j); %Nu tilde from previous iteration

    %Volumetric fluxes at faces 
    fW=-lgt_fw*u_face(i,j);
    fN=-lgt_fn*v_face(i,j); 
    fE= lgt_fe*u_face(i,j+1); 
    fS= lgt_fs*v_face(i+1,j);

    %CROSS DIFFUSION EXPLICIT COMPUTATION
    %Get nu_tilde vertexes values
    
    %Split nu_tildes
    nu_tilde_wn=nu_tilde_corners(i,j,1,1);
    nu_tilde_en=nu_tilde_corners(i,j,1,2);
    nu_tilde_es=nu_tilde_corners(i,j,1,3);
    nu_tilde_ws=nu_tilde_corners(i,j,1,4);

    %compute cross diffusion terms for each face 
    Sdnu_t_w=dW_c(i,j)*(nu_tilde_wn-nu_tilde_ws);%face w
    Sdnu_t_n=dN_c(i,j)*(nu_tilde_en-nu_tilde_wn);%face n
    Sdnu_t_e=dE_c(i,j)*(nu_tilde_en-nu_tilde_es);%face e
    Sdnu_t_s=dS_c(i,j)*(nu_tilde_es-nu_tilde_ws);%face s

    %Sum all contributions
    Sdnu_t_T= Sdnu_t_w+ Sdnu_t_n + Sdnu_t_e + Sdnu_t_s;

    %Nonlinear Diffusion 
    dnu_t_dx=grad_nu_tilde(i,j,1,1);
    dnu_t_dy=grad_nu_tilde(i,j,1,2);
  
    nonLinDiff=(cb2/sigma_sa)*(dnu_t_dx^2 + dnu_t_dy^2);

    %LINEARIZATION OF SOURCE TERMS FOR IMPLICIT APROXIMATION
    %  Implicit Taylor 1st order expansion
    %Q(nu_tilde_k+1)=~ Q(nu_tilde_k) + Q'(nu_tilde_k)*(nu_tilde_k+1
    % - nu_tilde_k)
    %central coeffitient
    %AP=-Q'(nu_tilde_k)*cellVol
    %Source
    %AS=(Q(nu_tilde_k) - Q'(nu_tilde_k)*nu_tilde_k))*cellVol

    %Q(nu_tilde)=Tp(nu_tilde) + Td(nu_tilde)

    d_wall=distMinWall(i,j);

    %__________ Turbulence Production Tp(nu_tilde) ___________

    %Viscous ratio*
    x_nu=nu_tilde_k/nu;
    x_nu_prime=1/nu;
    
    %%Viscous damping function 1 fv1*
    fv1=(x_nu^3)/(x_nu^3 + cv1_3);
    fv1_prime=(3*x_nu_prime*(x_nu^2)*cv1_3)/((x_nu^3 + cv1_3)^2);

    %Viscous damping function 2 fv2*
    fv2=1 - x_nu/(1 + x_nu*fv1);
    fv2_prime=(((x_nu^2)*fv1_prime) - x_nu_prime)/((1 + x_nu*fv1)^2);
     
    %abs  vorticity (gamma tilde)*
    s_vort=lower_omega(i,j);
    
    %modified vorticity*           
    s_vort_tilde = s_vort + (nu_tilde_k*fv2)/((kappa*d_wall)^2);
    s_vort_tilde_prime= (1/((kappa*d_wall)^2))*(fv2 +...
        fv2_prime*nu_tilde_k);
    
    %tp Turbulence production*
    tp=cb1*s_vort_tilde*nu_tilde_k;
    tp_prime= cb1*(s_vort_tilde_prime*nu_tilde_k + s_vort_tilde);
    
    %__________ Turbulence Destruction Td(nu_tilde) __________
    
    %function r*
    r=min(10,nu_tilde_k/(s_vort_tilde * (kappa*d_wall)^2 ));
    r_prime=(1/((s_vort_tilde*kappa*d_wall)^2))*(s_vort_tilde -...
        nu_tilde_k*s_vort_tilde_prime);
    
    %function g*
    g=r + cw2*(r^6 -r);
    g_prime=r_prime*(1-cw2 + 6*cw2*r^5);
    
    %wall damping function fw
    fw=g*(((1 + cw3_6)/(g^6 + cw3_6))^(1/6));
    fw_prime=g_prime*cw3_6*(((1 + cw3_6)/((g^6 + cw3_6)^7))^(1/6));

    %Td
    td=-cw1*fw*(nu_tilde_k/d_wall)^2;
    td_prime=-cw1*(fw_prime*(nu_tilde_k/d_wall)^2 + ...
        2*fw*(nu_tilde_k/(d_wall^2)));
    
    %ALL SOURCES Q(nu_tilde_k) & Q'(nu_tilde_k)
    q_s=tp + td;
    q_s_prime=tp_prime + td_prime;

    %Coeffitiens using upwind squeme

    Ant_W(i,j)=dW(i,j) + max(0,-fW);
    Ant_N(i,j)=dN(i,j) + max(0,-fE);
    Ant_E(i,j)=0;
    Ant_S(i,j)=dS(i,j);

    Ant_P(i,j)=Ant_W(i,j) + Ant_N(i,j)  +dS(i,j) +fW +fN +fE  - ...
        cellVols(i,j)*q_s_prime;

    Su_nt(i,j)=Sdnu_t_T + cellVols(i,j)*(nonLinDiff + q_s - ...
        q_s_prime*nu_tilde_k) + 0.1*nu*(dS(i,j)- fS) ;

    %west south_________________________________________________________
    i=ny;
    j=1;

    %DIFFUSION AND CONVECTION
    %get face areas
    lgt_fw=lgtFaces(i,j,1);
    lgt_fn=lgtFaces(i,j,2); 
    lgt_fe=lgtFaces(i,j,3); 
    lgt_fs=lgtFaces(i,j,4);

    nu_tilde_k=nu_tilde(i,j); %Nu tilde from previous iteration

    %Volumetric fluxes at faces 
    fW=-lgt_fw*u_face(i,j);
    fN=-lgt_fn*v_face(i,j); 
    fE= lgt_fe*u_face(i,j+1);   
    fS= lgt_fs*v_face(i+1,j);

    %CROSS DIFFUSION EXPLICIT COMPUTATION
    %Get nu_tilde vertexes values
    
    %Split nu_tildes
    nu_tilde_wn=nu_tilde_corners(i,j,1,1);
    nu_tilde_en=nu_tilde_corners(i,j,1,2);
    nu_tilde_es=nu_tilde_corners(i,j,1,3);
    nu_tilde_ws=nu_tilde_corners(i,j,1,4);

    %compute cross diffusion terms for each face 
    Sdnu_t_w=dW_c(i,j)*(nu_tilde_wn-nu_tilde_ws);%face w
    Sdnu_t_n=dN_c(i,j)*(nu_tilde_en-nu_tilde_wn);%face n
    Sdnu_t_e=dE_c(i,j)*(nu_tilde_en-nu_tilde_es);%face e
    Sdnu_t_s=dS_c(i,j)*(nu_tilde_es-nu_tilde_ws);%face s

    %Sum all contributions
    Sdnu_t_T= Sdnu_t_w+ Sdnu_t_n + Sdnu_t_e + Sdnu_t_s;

    %Nonlinear Diffusion 
    dnu_t_dx=grad_nu_tilde(i,j,1,1);
    dnu_t_dy=grad_nu_tilde(i,j,1,2);
  
    nonLinDiff=(cb2/sigma_sa)*(dnu_t_dx^2 + dnu_t_dy^2);

    %LINEARIZATION OF SOURCE TERMS FOR IMPLICIT APROXIMATION
    %  Implicit Taylor 1st order expansion
    %Q(nu_tilde_k+1)=~ Q(nu_tilde_k) + Q'(nu_tilde_k)*(nu_tilde_k+1
    % - nu_tilde_k)
    %central coeffitient
    %AP=-Q'(nu_tilde_k)*cellVol
    %Source
    %AS=(Q(nu_tilde_k) - Q'(nu_tilde_k)*nu_tilde_k))*cellVol

    %Q(nu_tilde)=Tp(nu_tilde) + Td(nu_tilde)

    d_wall=distMinWall(i,j);

    %__________ Turbulence Production Tp(nu_tilde) ___________

    %Viscous ratio*
    x_nu=nu_tilde_k/nu;
    x_nu_prime=1/nu;
    
    %%Viscous damping function 1 fv1*
    fv1=(x_nu^3)/(x_nu^3 + cv1_3);
    fv1_prime=(3*x_nu_prime*(x_nu^2)*cv1_3)/((x_nu^3 + cv1_3)^2);

    %Viscous damping function 2 fv2*
    fv2=1 - x_nu/(1 + x_nu*fv1);
    fv2_prime=(((x_nu^2)*fv1_prime) - x_nu_prime)/((1 + x_nu*fv1)^2);
     
    %abs  vorticity (gamma tilde)*
    s_vort=lower_omega(i,j);
    
    %modified vorticity*           
    s_vort_tilde = s_vort + (nu_tilde_k*fv2)/((kappa*d_wall)^2);
    s_vort_tilde_prime= (1/((kappa*d_wall)^2))*(fv2 +...
        fv2_prime*nu_tilde_k);
    
    %tp Turbulence production*
    tp=cb1*s_vort_tilde*nu_tilde_k;
    tp_prime= cb1*(s_vort_tilde_prime*nu_tilde_k + s_vort_tilde);
    
    %__________ Turbulence Destruction Td(nu_tilde) __________
    
    %function r*
    r=min(10,nu_tilde_k/(s_vort_tilde * (kappa*d_wall)^2 ));
    r_prime=(1/((s_vort_tilde*kappa*d_wall)^2))*(s_vort_tilde -...
        nu_tilde_k*s_vort_tilde_prime);
    
    %function g*
    g=r + cw2*(r^6 -r);
    g_prime=r_prime*(1-cw2 + 6*cw2*r^5);
    
    %wall damping function fw
    fw=g*(((1 + cw3_6)/(g^6 + cw3_6))^(1/6));
    fw_prime=g_prime*cw3_6*(((1 + cw3_6)/((g^6 + cw3_6)^7))^(1/6));

    %Td
    td=-cw1*fw*(nu_tilde_k/d_wall)^2;
    td_prime=-cw1*(fw_prime*(nu_tilde_k/d_wall)^2 + ...
        2*fw*(nu_tilde_k/(d_wall^2)));
    
    %ALL SOURCES Q(nu_tilde_k) & Q'(nu_tilde_k)
    q_s=tp + td;
    q_s_prime=tp_prime + td_prime;

    %Coeffitiens using upwind squeme
    Ant_W(i,j)=0;
    Ant_N(i,j)=dN(i,j) + max(0,-fE);
    Ant_E(i,j)=dE(i,j) + max(0,-fN);
    Ant_S(i,j)=0;

    Ant_P(i,j)=dW(i,j)+ Ant_N(i,j) + Ant_E(i,j) +dS(i,j) ...
        +fN +fE - cellVols(i,j)*q_s_prime;

    Su_nt(i,j)=Sdnu_t_T + cellVols(i,j)*(nonLinDiff + q_s - ...
        q_s_prime*nu_tilde_k) +0.1*nu*(dW(i,j) +dS(i,j) -fW -fS);

end

