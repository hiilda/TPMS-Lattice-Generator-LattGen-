function [out_p,x,y,z] = stlinput(inp,gridpoints)
st = stlread(inp);
% get stl vertices
stp = st.Points; stf = st.ConnectivityList;
% get bounding box
stp_xmax = max(stp(:,1)); stp_xmin = min(stp(:,1));
stp_ymax = max(stp(:,2)); stp_ymin = min(stp(:,2));
stp_zmax = max(stp(:,3)); stp_zmin = min(stp(:,3));

%generate gridpoints for x, y, and z
x = linspace(stp_xmin, stp_xmax, gridpoints); 
y = linspace(stp_ymin, stp_ymax, gridpoints); 
z = linspace(stp_zmin, stp_zmax, gridpoints);

%create meshgrid for x, y, and z
[x1,y1,z1] = meshgrid(x,y,z);

out = intriangulation(stp,stf,[x1(:),y1(:),z1(:)]); out(out<0) = 0;
out_p = reshape(out,gridpoints,gridpoints,gridpoints);




