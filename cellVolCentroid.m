function [cell_centroid_x,cell_centroid_y,volcel] = ...
    cellVolCentroid(x1,y1,x2,y2,x3,y3,x4,y4)
%Compute the volume of the cell and the centroid dividing the cell in two
%triangles 

%cell area and centroid
%vertexes of the cell
p1=[x1,y1,0];
p2=[x2,y2,0];
p3=[x3,y3,0];
p4=[x4,y4,0];

%Dividing the cell in 2 triangles by creating 3 vectors
v1=p2-p1;
v2=p3-p1;
v3=p4-p1;

%Area of triangle 1
Area_1=0.5*vecnorm(cross(v1,v2));

%Area of triangle 2
Area_2=0.5*vecnorm(cross(v2,v3));

%Area of the cell
volcel =  Area_1 + Area_2;

%centroid of triangle 1
t1_cent_x=(p1(1) + p2(1) +p3(1))/3;
t1_cent_y=(p1(2) + p2(2) +p3(2))/3;

%Centroid of triangle 2
t2_cent_x=(p1(1) + p3(1) +p4(1))/3;
t2_cent_y=(p1(2) + p3(2) +p4(2))/3;

%centroid of the cell
cell_centroid_x=(t1_cent_x*Area_1 + t2_cent_x*Area_2)/volcel;
cell_centroid_y=(t1_cent_y*Area_1 + t2_cent_y*Area_2)/volcel;

end