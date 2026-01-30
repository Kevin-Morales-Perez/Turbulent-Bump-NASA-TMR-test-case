function [rsid_cont,err_cont] = checkContinuity(u_f,v_f,lgtFaces,nx,ny)
    %CONTINUITY CHECK

    rsid_cont=zeros(ny,nx);

    for i=1:ny
        for j=1:nx
            %Get cell faces areas
            lgt_fw=lgtFaces(i,j,1);%West
            lgt_fn=lgtFaces(i,j,2);%North 
            lgt_fe=lgtFaces(i,j,3); %East
            lgt_fs=lgtFaces(i,j,4);%South

            rsid_cont(i,j)=(lgt_fe*u_f(i,j+1) -lgt_fw*u_f(i,j))  + ...
                (-lgt_fn*v_f(i,j) + lgt_fs*v_f(i+1,j));

        end
    end

    err_cont=rms(rsid_cont(:));


end