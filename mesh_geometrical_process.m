function [cellVertxs,cellCentrs,lgtFaces,uVecParlFaces,faceCentrs,...
    cellVols,uVecNormFaces,distNeighbNods,uVecNeighbNods,...
    wlsqOperator,VecCentNodFaceCents,VecCentVertx,VecCentNbNods,...
    weightDistFactors,Xctrs,Yctrs,distMinWall] =...
    mesh_geometrical_process(X,Y,nx,ny,cellVertxs,cellCentrs,...
    lgtFaces,uVecParlFaces,faceCentrs,cellVols,uVecNormFaces,...
    distNeighbNods,uVecNeighbNods,wlsqOperator,VecCentNodFaceCents,...
    VecCentVertx,VecCentNbNods,weightDistFactors,Xctrs,Yctrs,distMinWall,...
    distPlat)

    %Compute cell vertex
    
    for i=1:ny %Iterate over rows
        for j=1:nx %Iterate over colums

            %Compute cell vertex
            x1=X(i,j);      y1=Y(i,j);
            x2=X(i,j+1);    y2=Y(i,j+1);
            x3=X(i+1,j+1);  y3=Y(i+1,j+1);
            x4=X(i+1,j);    y4=Y(i+1,j);

            Vwn=[x1,y1];%vertex west north
            Ven=[x2,y2];%vertex east north
            Ves=[x3,y3];%vertex east south
            Vws=[x4,y4];%vertex west south

            cellVertxs(i,j,:,:)=[Vwn;Ven;Ves;Vws];

            %_____________ New procedure for cell centers _________________
            
            [xc,yc,celvol]=cellVolCentroid(x1,y1,x2,y2,x3,y3,x4,y4);

            cell_centroid=[xc,yc];
            
            cellCentrs(i,j,:)=cell_centroid;

            %___________________________________________________________

            %Use the cell vertexs to create vectors parallel to faces
            %VECTORS FOR FACE W AND E ARE ALLWAYS POSITIVE AND
            %ONLY WITH Y COMPONEN
            %VECTORS FOR  FACES N AND S HAVE ALWAYS POSITIVE X COMPONENT

            par_w=Vwn-Vws;
            par_n=Ven-Vwn;
            par_e=Ven-Ves;
            par_s=Ves-Vws;

            %calculate the length of each face

            lgt_w=vecnorm(par_w);
            lgt_n=vecnorm(par_n);
            lgt_e=vecnorm(par_e);
            lgt_s=vecnorm(par_s);
            %Store variables
            lgtFaces(i,j,:) = [lgt_w, lgt_n, lgt_e, lgt_s];

            %make unitary previous vectors;
            upar_w=par_w/lgt_w;
            upar_n=par_n/lgt_n;
            upar_e=par_e/lgt_e;
            upar_s=par_s/lgt_s;

            %Store unitary vectors paralel to faces
            uVecParlFaces(i,j,:,:)=[upar_w;upar_n;upar_e;upar_s];

            %Center of faces 

            %face w, starting point:vertex vws, 
            
            fcent_w=Vws +0.5*lgt_w*upar_w;
            fcent_n=Vwn +0.5*lgt_n*upar_n;
            fcent_e=Ves +0.5*lgt_e*upar_e;
            fcent_s=Vws +0.5*lgt_s*upar_s;

            faceCentrs(i,j,:,:)=[fcent_w;fcent_n;fcent_e;fcent_s];

            %Cell Vols
            cellVols(i,j)=celvol;

            %___________________________________
            %unitary vectors normal to faces
            unorm_w=[-upar_w(2),upar_w(1)];
            unorm_n=[-upar_n(2),upar_n(1)];
            unorm_e=[upar_e(2),-upar_e(1)];
            unorm_s=[upar_s(2),-upar_s(1)];

            uVecNormFaces(i,j,:,:)=[unorm_w;unorm_n;unorm_e;unorm_s];

            %Vector from central node to center of faces

            %caculate vectors
            vec_p_fw=fcent_w-cell_centroid;
            vec_p_fn=fcent_n-cell_centroid;
            vec_p_fe=fcent_e-cell_centroid;
            vec_p_fs=fcent_s-cell_centroid;
            

            %Store vectors
            VecCentNodFaceCents(i,j,:,:)=[vec_p_fw;vec_p_fn;vec_p_fe;...
                vec_p_fs];


        end
    end

    %Variables that require previously computed values 


    %INTERIOR CELLS 
    %distance between central node and neigborhood nodes
    %unitary vectors from central node to neigborhood nodes
   

    for i=2:ny-1 %Iterate over rows
        for j=2:nx-1 %Iterate over columns

            %get the centroids of the nodes 
            cent_p=reshape(cellCentrs(i,j,:),[1,2]);
            cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
            cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
            cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
            cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);

            vec_p_w=cent_w-cent_p;
            vec_p_n=cent_n-cent_p;
            vec_p_e=cent_e-cent_p;
            vec_p_s=cent_s-cent_p;

            %Get the distances between nodes 
            lgt_p_w=vecnorm(vec_p_w);
            lgt_p_n=vecnorm(vec_p_n);
            lgt_p_e=vecnorm(vec_p_e);
            lgt_p_s=vecnorm(vec_p_s);

            %STORE the  distances
            distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];


            %Calculate unitary vectors 
            uvec_p_w=vec_p_w/lgt_p_w;
            uvec_p_n=vec_p_n/lgt_p_n;
            uvec_p_e=vec_p_e/lgt_p_e;
            uvec_p_s=vec_p_s/lgt_p_s;

            %STORE the values 

            uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];

            %WEIGHTED LEAST SQUARE OPERATOR

            %Neighborhood nodes grouped in a vector
            ngb_nodes=[cent_w;cent_n;cent_e;cent_s];

            dist_pn=ngb_nodes - cent_p;
            %Vector of weights
            norm=vecnorm(dist_pn,2,2);
            w_v=1./norm;
            w_v=diag(w_v);
            %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
            G_m=dist_pn'*(w_v*dist_pn);
            gradient_weights_op=G_m\(dist_pn'*w_v);

            wlsqOperator(i,j,:,:)=gradient_weights_op;
         
        end
    end

   
    %If there is no neighborhood node is calculated
    %the distance to the center of the face instead

    %Sites where no neighboorhood nodes available
    
    %left  wall (west)
    j=1;
    for i =2:ny-1
        cent_p=reshape(cellCentrs(i,j,:),[1,2]);
        cent_w=reshape(faceCentrs(i,j,1,:),[1,2]);%face center
        cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
        cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
        cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);

        vec_p_w=cent_w-cent_p;
        vec_p_n=cent_n-cent_p;
        vec_p_e=cent_e-cent_p;
        vec_p_s=cent_s-cent_p;

        %Get the distances between nodes 
        lgt_p_w=vecnorm(vec_p_w);
        lgt_p_n=vecnorm(vec_p_n);
        lgt_p_e=vecnorm(vec_p_e);
        lgt_p_s=vecnorm(vec_p_s);

        %STORE the  distances
        distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];


        %Calculate unitary vectors 
        uvec_p_w=vec_p_w/lgt_p_w;
        uvec_p_n=vec_p_n/lgt_p_n;
        uvec_p_e=vec_p_e/lgt_p_e;
        uvec_p_s=vec_p_s/lgt_p_s;

        %STORE the values 

        uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];
       
    end

    %Top wall north
    i=1;
    
    for j=2:nx-1
        cent_p=reshape(cellCentrs(i,j,:),[1,2]);
        cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
        cent_n=reshape(faceCentrs(i,j,2,:),[1,2]);%face center
        cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
        cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);

        vec_p_w=cent_w-cent_p;
        vec_p_n=cent_n-cent_p;
        vec_p_e=cent_e-cent_p;
        vec_p_s=cent_s-cent_p;

        %Get the distances between nodes 
        lgt_p_w=vecnorm(vec_p_w);
        lgt_p_n=vecnorm(vec_p_n);
        lgt_p_e=vecnorm(vec_p_e);
        lgt_p_s=vecnorm(vec_p_s);

        %STORE the  distances
        distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];

        %Calculate unitary vectors 
        uvec_p_w=vec_p_w/lgt_p_w;
        uvec_p_n=vec_p_n/lgt_p_n;
        uvec_p_e=vec_p_e/lgt_p_e;
        uvec_p_s=vec_p_s/lgt_p_s;

        %STORE the values 

        uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];

    end

    %1.2.3- Right Wall - East 
    
    j=nx;
    for i=2:ny-1

        cent_p=reshape(cellCentrs(i,j,:),[1,2]);
        cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
        cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
        cent_e=reshape(faceCentrs(i,j,3,:),[1,2]);%face center
        cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);

        vec_p_w=cent_w-cent_p;
        vec_p_n=cent_n-cent_p;
        vec_p_e=cent_e-cent_p;
        vec_p_s=cent_s-cent_p;

        %Get the distances between nodes 
        lgt_p_w=vecnorm(vec_p_w);
        lgt_p_n=vecnorm(vec_p_n);
        lgt_p_e=vecnorm(vec_p_e);
        lgt_p_s=vecnorm(vec_p_s);

        %STORE the  distances
        distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];


        %Calculate unitary vectors 
        uvec_p_w=vec_p_w/lgt_p_w;
        uvec_p_n=vec_p_n/lgt_p_n;
        uvec_p_e=vec_p_e/lgt_p_e;
        uvec_p_s=vec_p_s/lgt_p_s;

        %STORE the values 

        uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];        

    end

    %1.2.4- Bottom wall -- South
    
    i=ny;
    
    for j=2:nx-1
        cent_p=reshape(cellCentrs(i,j,:),[1,2]);
        cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
        cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
        cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
        cent_s=reshape(faceCentrs(i,j,4,:),[1,2]);%face center

        vec_p_w=cent_w-cent_p;
        vec_p_n=cent_n-cent_p;
        vec_p_e=cent_e-cent_p;
        vec_p_s=cent_s-cent_p;

        %Get the distances between nodes 
        lgt_p_w=vecnorm(vec_p_w);
        lgt_p_n=vecnorm(vec_p_n);
        lgt_p_e=vecnorm(vec_p_e);
        lgt_p_s=vecnorm(vec_p_s);

        %STORE the  distances
        distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];


        %Calculate unitary vectors 
        uvec_p_w=vec_p_w/lgt_p_w;
        uvec_p_n=vec_p_n/lgt_p_n;
        uvec_p_e=vec_p_e/lgt_p_e;
        uvec_p_s=vec_p_s/lgt_p_s;

        %STORE the values 

        uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];

    end

    %CORNERS
    %1.3.1- North-West
    i=1;
    j=1;

    cent_p=reshape(cellCentrs(i,j,:),[1,2]);
    cent_w=reshape(faceCentrs(i,j,1,:),[1,2]);%face center
    cent_n=reshape(faceCentrs(i,j,2,:),[1,2]);%face center
    cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
    cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);

    vec_p_w=cent_w-cent_p;
    vec_p_n=cent_n-cent_p;
    vec_p_e=cent_e-cent_p;
    vec_p_s=cent_s-cent_p;

    %Get the distances between nodes 
    lgt_p_w=vecnorm(vec_p_w);
    lgt_p_n=vecnorm(vec_p_n);
    lgt_p_e=vecnorm(vec_p_e);
    lgt_p_s=vecnorm(vec_p_s);

    %STORE the  distances
    distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];

    %Calculate unitary vectors 
    uvec_p_w=vec_p_w/lgt_p_w;
    uvec_p_n=vec_p_n/lgt_p_n;
    uvec_p_e=vec_p_e/lgt_p_e;
    uvec_p_s=vec_p_s/lgt_p_s;

    %STORE the values 

    uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];

    %1.3.2- North-East
    i=1;
    j=nx;

    cent_p=reshape(cellCentrs(i,j,:),[1,2]);
    cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
    cent_n=reshape(faceCentrs(i,j,2,:),[1,2]);%face center
    cent_e=reshape(faceCentrs(i,j,3,:),[1,2]);%face center
    cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);

    vec_p_w=cent_w-cent_p;
    vec_p_n=cent_n-cent_p;
    vec_p_e=cent_e-cent_p;
    vec_p_s=cent_s-cent_p;

    %Get the distances between nodes 
    lgt_p_w=vecnorm(vec_p_w);
    lgt_p_n=vecnorm(vec_p_n);
    lgt_p_e=vecnorm(vec_p_e);
    lgt_p_s=vecnorm(vec_p_s);

    %STORE the  distances
    distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];


    %Calculate unitary vectors 
    uvec_p_w=vec_p_w/lgt_p_w;
    uvec_p_n=vec_p_n/lgt_p_n;
    uvec_p_e=vec_p_e/lgt_p_e;
    uvec_p_s=vec_p_s/lgt_p_s;

    %STORE the values 

    uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];

    %1.3.3- South-East (Symetric at South and outlet at East)
    i=ny;
    j=nx;

    cent_p=reshape(cellCentrs(i,j,:),[1,2]);
    cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
    cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
    cent_e=reshape(faceCentrs(i,j,3,:),[1,2]);%face center
    cent_s=reshape(faceCentrs(i,j,4,:),[1,2]);%face center

    vec_p_w=cent_w-cent_p;
    vec_p_n=cent_n-cent_p;
    vec_p_e=cent_e-cent_p;
    vec_p_s=cent_s-cent_p;

    %Get the distances between nodes 
    lgt_p_w=vecnorm(vec_p_w);
    lgt_p_n=vecnorm(vec_p_n);
    lgt_p_e=vecnorm(vec_p_e);
    lgt_p_s=vecnorm(vec_p_s);

    %STORE the  distances
    distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];

    %Calculate unitary vectors 
    uvec_p_w=vec_p_w/lgt_p_w;
    uvec_p_n=vec_p_n/lgt_p_n;
    uvec_p_e=vec_p_e/lgt_p_e;
    uvec_p_s=vec_p_s/lgt_p_s;

    %STORE the values 

    uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];

    %1.3.4- South-West
    i=ny;
    j=1;

    cent_p=reshape(cellCentrs(i,j,:),[1,2]);
    cent_w=reshape(faceCentrs(i,j,1,:),[1,2]);%face center
    cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
    cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
    cent_s=reshape(faceCentrs(i,j,4,:),[1,2]);%face center

    vec_p_w=cent_w-cent_p;
    vec_p_n=cent_n-cent_p;
    vec_p_e=cent_e-cent_p;
    vec_p_s=cent_s-cent_p;

    %Get the distances between nodes 
    lgt_p_w=vecnorm(vec_p_w);
    lgt_p_n=vecnorm(vec_p_n);
    lgt_p_e=vecnorm(vec_p_e);
    lgt_p_s=vecnorm(vec_p_s);

    %STORE the  distances
    distNeighbNods(i,j,:)=[lgt_p_w,lgt_p_n,lgt_p_e,lgt_p_s];

    %Calculate unitary vectors 
    uvec_p_w=vec_p_w/lgt_p_w;
    uvec_p_n=vec_p_n/lgt_p_n;
    uvec_p_e=vec_p_e/lgt_p_e;
    uvec_p_s=vec_p_s/lgt_p_s;

    %STORE the values 

    uVecNeighbNods(i,j,:,:)=[uvec_p_w;uvec_p_n;uvec_p_e;uvec_p_s];

    %______________________________________________________________
    %Additional points in faces for weighted least square gradients
    %West Edge 

    %aditional centers for least square gradient west edge
    centPointAddBoundWest=zeros(ny,2);
    %aditional centers for least square gradient North edge
    centPointAddBoundNorth=zeros(nx,2);
    %aditional centers for least square gradient East edge
    centPointAddBoundEast=zeros(ny,2);
    %aditional centers for least square gradient south edge
    centPointAddBoundSouth=zeros(nx,2);


    j=1;
    for i=1:ny
        
        centPointAddBoundWest(i,:)=reshape(faceCentrs(i,j,1,:),[1,2]);
    end
    %North
    i=1;
    for j=1:nx
        centPointAddBoundNorth(j,:)=reshape(faceCentrs(i,j,2,:),[1,2]);
    end
    %East
    j=nx;
    for i=1:ny
        centPointAddBoundEast(i,:)=reshape(faceCentrs(i,j,3,:),[1,2]);
    end
    %South (Special treatment to simulate a perpendicular line to the
    %solid surface
    i=ny;
    for j =1:nx

        %Get vertex west south from cell
        vert_ws= reshape(cellVertxs(i,j,4,:),[1,2]);
        %get cell center
        cel_cent=reshape(cellCentrs(i,j,:),[1,2]);
        %get unitary vector tangential to face S
        eta_s =reshape(uVecParlFaces(i,j,4,:),[1,2]);

        %calculate point coordinates for wall pressure boundary
        %the line that joints this point to cell center is perpendicular 
        %to face a, then least weighted squares gradient can be used to 
        % compute
        %pressure gradient at wall with non orthogonal faces
        centPointAddBoundSouth(j,:)=vert_ws +...
            dot((cel_cent-vert_ws),eta_s)*eta_s;

    end

    %Replace face center at south for centPointAddBoundSouth
    i=ny;
    for j=1:nx
        faceCentrs(i,j,4,:)=centPointAddBoundSouth(j,:);
    end

    %Additional Weighted Least Square Gradients

    %west edge
    j=1;
    for i = 2:ny-1
        %get the centroids of the nodes (replace west node to additional
        %center west face node
        cent_p=reshape(cellCentrs(i,j,:),[1,2]);
        cent_w=centPointAddBoundWest(i,:);
        cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
        cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
        cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);

        %WEIGHTED LEAST SQUARE OPERATOR

        %Neighborhood nodes grouped in a vector
        ngb_nodes=[cent_w;cent_n;cent_e;cent_s];

        dist_pn=ngb_nodes - cent_p;
        %Vector of weights
        norm=vecnorm(dist_pn,2,2);
        w_v=1./norm;
        w_v=diag(w_v);
        %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
        G_m=dist_pn'*(w_v*dist_pn);
        gradient_weights_op=G_m\(dist_pn'*w_v);

        wlsqOperator(i,j,:,:)=gradient_weights_op;

    end
    
    %North Edge
    i=1;
    for j=2:nx-1
        %node center coordinates
        cent_p=reshape(cellCentrs(i,j,:),[1,2]);
        cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
        cent_n=centPointAddBoundNorth(j,:);
        cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
        cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);
        
        %WEIGHTED LEAST SQUARE OPERATOR
        
        %Neighborhood nodes grouped in a vector
        ngb_nodes=[cent_w;cent_n;cent_e;cent_s];
        
        dist_pn=ngb_nodes - cent_p;
        %Vector of weights
        norm=vecnorm(dist_pn,2,2);
        w_v=1./norm;
        w_v=diag(w_v);
        %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
        G_m=dist_pn'*(w_v*dist_pn);
        gradient_weights_op=G_m\(dist_pn'*w_v);
        
        wlsqOperator(i,j,:,:)=gradient_weights_op;

    end
    
    %East Edge 
    j=nx;
    for i =2:ny-1
        %node center coordinates
        cent_p=reshape(cellCentrs(i,j,:),[1,2]);
        cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
        cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
        cent_e=centPointAddBoundEast(i,:);
        cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);
        
        %WEIGHTED LEAST SQUARE OPERATOR
        
        %Neighborhood nodes grouped in a vector
        ngb_nodes=[cent_w;cent_n;cent_e;cent_s];
        
        dist_pn=ngb_nodes - cent_p;
        %Vector of weights
        norm=vecnorm(dist_pn,2,2);
        w_v=1./norm;
        w_v=diag(w_v);
        %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
        G_m=dist_pn'*(w_v*dist_pn);
        gradient_weights_op=G_m\(dist_pn'*w_v);
        
        wlsqOperator(i,j,:,:)=gradient_weights_op;


    end

    %South Edge 

    i=ny;
    for j=2:nx-1
        %node center coordinates
        cent_p=reshape(cellCentrs(i,j,:),[1,2]);
        cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
        cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
        cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
        cent_s=centPointAddBoundSouth(j,:);
        
        %WEIGHTED LEAST SQUARE OPERATOR
        
        %Neighborhood nodes grouped in a vector
        ngb_nodes=[cent_w;cent_n;cent_e;cent_s];
        
        dist_pn=ngb_nodes - cent_p;
        %Vector of weights
        norm=vecnorm(dist_pn,2,2);
        w_v=1./norm;
        w_v=diag(w_v);
        %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
        G_m=dist_pn'*(w_v*dist_pn);
        gradient_weights_op=G_m\(dist_pn'*w_v);
        
        wlsqOperator(i,j,:,:)=gradient_weights_op;

    end

    %West north corner
    i=1;
    j=1;

    %node center coordinates
    cent_p=reshape(cellCentrs(i,j,:),[1,2]);
    cent_w=centPointAddBoundWest(i,:);
    cent_n=centPointAddBoundNorth(j,:);
    cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
    cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);
    
    %WEIGHTED LEAST SQUARE OPERATOR
    
    %Neighborhood nodes grouped in a vector
    ngb_nodes=[cent_w;cent_n;cent_e;cent_s];
    
    dist_pn=ngb_nodes - cent_p;
    %Vector of weights
    norm=vecnorm(dist_pn,2,2);
    w_v=1./norm;
    w_v=diag(w_v);
    %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
    G_m=dist_pn'*(w_v*dist_pn);
    gradient_weights_op=G_m\(dist_pn'*w_v);
    
    wlsqOperator(i,j,:,:)=gradient_weights_op;



    %east north corner
    i=1;
    j=nx;

    %node center coordinates
    cent_p=reshape(cellCentrs(i,j,:),[1,2]);
    cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
    cent_n=centPointAddBoundNorth(j,:);
    cent_e=centPointAddBoundEast(i,:);
    cent_s=reshape(cellCentrs(i+1,j,:),[1,2]);
    
    %WEIGHTED LEAST SQUARE OPERATOR
    
    %Neighborhood nodes grouped in a vector
    ngb_nodes=[cent_w;cent_n;cent_e;cent_s];
    
    dist_pn=ngb_nodes - cent_p;
    %Vector of weights
    norm=vecnorm(dist_pn,2,2);
    w_v=1./norm;
    w_v=diag(w_v);
    %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
    G_m=dist_pn'*(w_v*dist_pn);
    gradient_weights_op=G_m\(dist_pn'*w_v);
    
    wlsqOperator(i,j,:,:)=gradient_weights_op;

    %east south corner
    i=ny;
    j=nx;

    %node center coordinates
    cent_p=reshape(cellCentrs(i,j,:),[1,2]);
    cent_w=reshape(cellCentrs(i,j-1,:),[1,2]);
    cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
    cent_e=centPointAddBoundEast(i,:);
    cent_s=centPointAddBoundSouth(j,:);
    
    %WEIGHTED LEAST SQUARE OPERATOR
    
    %Neighborhood nodes grouped in a vector
    ngb_nodes=[cent_w;cent_n;cent_e;cent_s];
    
    dist_pn=ngb_nodes - cent_p;
    %Vector of weights
    norm=vecnorm(dist_pn,2,2);
    w_v=1./norm;
    w_v=diag(w_v);
    %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
    G_m=dist_pn'*(w_v*dist_pn);
    gradient_weights_op=G_m\(dist_pn'*w_v);
    
    wlsqOperator(i,j,:,:)=gradient_weights_op;

    %west south corner
    i=ny;
    j=1;

    %node center coordinates
    cent_p=reshape(cellCentrs(i,j,:),[1,2]);
    cent_w=centPointAddBoundWest(i,:);
    cent_n=reshape(cellCentrs(i-1,j,:),[1,2]);
    cent_e=reshape(cellCentrs(i,j+1,:),[1,2]);
    cent_s=centPointAddBoundSouth(j,:);
    
    %WEIGHTED LEAST SQUARE OPERATOR
    
    %Neighborhood nodes grouped in a vector
    ngb_nodes=[cent_w;cent_n;cent_e;cent_s];
    
    dist_pn=ngb_nodes - cent_p;
    %Vector of weights
    norm=vecnorm(dist_pn,2,2);
    w_v=1./norm;
    w_v=diag(w_v);
    %g matrix , this are the matrices that multiplies the gradient at the left side of the equation
    G_m=dist_pn'*(w_v*dist_pn);
    gradient_weights_op=G_m\(dist_pn'*w_v);
    
    wlsqOperator(i,j,:,:)=gradient_weights_op;

    %_______________________vertexes vector_______________________________

    %vectors from cell center each to vertex
     for i=1:ny
        for j=1:nx

            %get cell vertex
            cell_vertex=reshape(cellVertxs(i,j,:,:),[4,2]);
            vrt_wn=cell_vertex(1,:);
            vrt_en=cell_vertex(2,:);
            vrt_es=cell_vertex(3,:);
            vrt_ws=cell_vertex(4,:);

            %get cell centroid
            cell_cent=reshape(cellCentrs(i,j,:),[1,2]);

            %compute vectors;
            vec_c_wn=vrt_wn-cell_cent;
            vec_c_en=vrt_en-cell_cent;
            vec_c_es=vrt_es-cell_cent;
            vec_c_ws=vrt_ws-cell_cent;

            %Store values 

            VecCentVertx(i,j,:,:)=[vec_c_wn;vec_c_en;vec_c_es;vec_c_ws];

        end
     end

     %_________________Weight factors____________________________________

     for i=1:ny
         for j=1:nx

             %get distances from cell center to neighborhood nodes
             dist_nods=reshape(distNeighbNods(i,j,:,:),[1,4]);

             %split dist_nods
             dist_p_w=dist_nods(1);
             dist_p_n=dist_nods(2);
             dist_p_e=dist_nods(3);
             dist_p_s=dist_nods(4);

             %get distances from cell center to each face center

             %get vectors from cell center to face centers

             vec_cent_faces=reshape(VecCentNodFaceCents(i,j,:,:),[4,2]);

             %split vec_cent_faces

             vec_cent_w=vec_cent_faces(1,:);
             vec_cent_n=vec_cent_faces(2,:);
             vec_cent_e=vec_cent_faces(3,:);
             vec_cent_s=vec_cent_faces(4,:);

             %get the lenght of the vectors

             dist_cent_fw=vecnorm(vec_cent_w);
             dist_cent_fn=vecnorm(vec_cent_n);
             dist_cent_fe=vecnorm(vec_cent_e);
             dist_cent_fs=vecnorm(vec_cent_s);

             %compute weight factors
             weight_c_w=(dist_p_w - dist_cent_fw)/dist_p_w;
             weight_c_n=(dist_p_n - dist_cent_fn)/dist_p_n;
             weight_c_e=(dist_p_e - dist_cent_fe)/dist_p_e;
             weight_c_s=(dist_p_s - dist_cent_fs)/dist_p_s;

             %Store values

             weightDistFactors(i,j,:)=[weight_c_w,weight_c_n,weight_c_e...
                 ,weight_c_s];

         end
     end

     %___________________Mesh with cell centers__________________________
    for i=1:ny
        for j=1:nx
            node_cord=reshape(cellCentrs(i,j,:),[1,2]);
            Xctrs(i,j)=node_cord(1);
            Yctrs(i,j)=node_cord(2);

        end
    end
    %___________MINIMUM DISTANCE TO THE WALL COMPUTATION_____________________
   
    %k_1 = pi/0.9;
    %k_2 = k_1*(distPlat + 0.3);
   
    bump_start = distPlat + 0.3;
    bump_end = distPlat + 1.2;
    plate_end = distPlat + 1.5;  % Assuming bumpLgt = 1.5
   
    for i = 1:ny
        for j = 1:nx
            % Get the x & y coordinate of the center of the cell
            x_cent = cellCentrs(i,j,1);
            y_cent = cellCentrs(i,j,2);
            
            % Case 1: Point above the bump region
            if x_cent > bump_start && x_cent < bump_end
                %if i > ny-12


                    %cells very close to the wall

                    %use dot product aproximation

                    %get south west vertex of the cell ny,j

                    vert_sw=reshape(cellVertxs(ny,j,4,:),[1,2]);

                    %get south face normal unitary vector from cell ny,j

                    norm_s= reshape(uVecNormFaces(ny,j,4,:),[1,2]);

                    %vector from vertex sw to cell centroid

                    vec_cent_sw =[x_cent,y_cent] - vert_sw;

                    %dot vector with inverse unitary normal 

                    distMinWall(i,j)= dot(vec_cent_sw,(-norm_s)); 
                
                
                       
            % Case 2: Point above the flat portion of the plate
            elseif (x_cent >= distPlat && x_cent < bump_start) || ...
                   (x_cent >= bump_end && x_cent < plate_end)
                
                distMinWall(i,j) = y_cent;
                    
            % Case 3: Cell centroid before or after the plate
            else
                if x_cent < distPlat
                    % Point before the plate - distance to leading edge
                    p_wall = [distPlat, 0];
                else
                    % Point after the plate - distance to trailing edge
                    p_wall = [plate_end, 0];
                end
                distMinWall(i,j) = sqrt((x_cent - p_wall(1))^2 + (y_cent - p_wall(2))^2);
            end
        end
    end

    %______________________Vectors from cell center to NB nodes__________

    for i=1:ny
        for j=1:nx
            for k=1:4
                u_vec_p_a=reshape(uVecNeighbNods(i,j,k,:),[1,2]);
                dist_p_a = distNeighbNods(i,j,k);
                VecCentNbNods(i,j,k,:)=u_vec_p_a*dist_p_a;
            end
        end
    end

end