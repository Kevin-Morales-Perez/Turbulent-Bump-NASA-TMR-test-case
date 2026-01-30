%V velocity interpolation test 

v_face_test=zeros(ny+1,nx);

for i=1:ny-1
    for j=1:nx

        

        %Velocity at cell P
        u_p=u(i,j);
        v_p=v(i,j);

        U_P=[u_p,v_p];

        %Velocity at cell S
        u_s=u(i+1,j);
        v_s=v(i+1,j);

        U_S=[u_s,v_s];

        %Get unitary vector normal to face S
        uvecn_s=reshape(uVecNormFaces(i,j,4,:),[1,2]);

        %weight distance factor at face S
        weight_S=weightDistFactors(i,j,4);

        %Averaged velocity vector at face S 
        U_s=weight_S*U_P + (1-weight_S)*U_S;

        v_face_test(i+1,j)=dot(U_s,uvecn_s);


    end
end

v_face_test(1,:)=v(1,:);