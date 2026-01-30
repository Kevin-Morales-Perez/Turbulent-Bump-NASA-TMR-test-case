%cell area and centroid

%vertexes of the cell
p1=[0,5,0];
p2=[12,4,0];
p3=[12,-4,0];
p4=[0,0,0];

%Dividing the cell in 2 triangles by creating 3 vectors
v1=p2-p1;
v2=p3-p1;
v3=p4-p1;

%Area of triangle 1
Area_1=0.5*vecnorm(cross(v1,v2));

%Area of triangle 2
Area_2=0.5*vecnorm(cross(v2,v3));

%Area of the cell
Area_total =  Area_1 + Area_2;

%centroid of triangle 1
t1_cent_x=(p1(1) + p2(1) +p3(1))/3;
t1_cent_y=(p1(2) + p2(2) +p3(2))/3;

%Centroid of triangle 2
t2_cent_x=(p1(1) + p3(1) +p4(1))/3;
t2_cent_y=(p1(2) + p3(2) +p4(2))/3;

%centroid of the cell
cell_centroid_x=(t1_cent_x*Area_1 + t2_cent_x*Area_2)/Area_total;
cell_centroid_y=(t1_cent_y*Area_1 + t2_cent_y*Area_2)/Area_total;


















