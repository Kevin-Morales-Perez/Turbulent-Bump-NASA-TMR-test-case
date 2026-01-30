%Test Unitary Vectors Validity

valid_uVecnormsFaces_W=zeros(ny,nx);
valid_uVecnormsFaces_N=zeros(ny,nx);
valid_uVecnormsFaces_E=zeros(ny,nx);
valid_uVecnormsFaces_S=zeros(ny,nx);


tol_len =1e-6;

tol_matrix_W=zeros(ny,nx);
tol_matrix_N=zeros(ny,nx);
tol_matrix_E=zeros(ny,nx);
tol_matrix_S=zeros(ny,nx);




for i=1:ny
    for j=1:nx
        
        %_W_
        vec_x = reshape(uVecNormFaces(i,j,1,:),[1,2]);
        vec_x_length = vecnorm(vec_x);

        tol_matrix_W(i,j)=abs(1- vec_x_length);

        if tol_matrix_W(i,j) > tol_len
            valid_uVecnormsFaces_W(i,j)=1;
        end

        %_N_
        vec_x = reshape(uVecNormFaces(i,j,2,:),[1,2]);
        vec_x_length = vecnorm(vec_x);

        tol_matrix_N(i,j)=abs(1- vec_x_length);

        if tol_matrix_N(i,j) > tol_len
            valid_uVecnormsFaces_N(i,j)=1;
        end

        %_E_
        vec_x = reshape(uVecNormFaces(i,j,3,:),[1,2]);
        vec_x_length = vecnorm(vec_x);

        tol_matrix_E(i,j)=abs(1- vec_x_length);

        if tol_matrix_E(i,j) > tol_len
            valid_uVecnormsFaces_E(i,j)=1;
        end

        %_S_
        vec_x = reshape(uVecNormFaces(i,j,4,:),[1,2]);
        vec_x_length = vecnorm(vec_x);

        tol_matrix_S(i,j)=abs(1- vec_x_length);

        if tol_matrix_S(i,j) > tol_len
            valid_uVecnormsFaces_S(i,j)=1;
        end

        
    end
end
