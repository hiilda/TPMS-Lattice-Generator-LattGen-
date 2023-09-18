clc
clear all

% cuboid domain dimensions
len = 40; bre = 40; wid = 20;

grade = 'No'; equation_type = 'Circular'; %user selects type of grading to use 
positionX = 0.5; positionY = 0.5; positionZ = 0.5; %reference position for grading
rr = 40; % radius of circle if needed by user
m1 = 1; m2 = 0; m3 = 0; % coefficients for variables in equations
rd1 = 0.1; rd2 = 1; app_t1 = 0.4; app_t2 = 2; % relative densities and approximate thicknesses for grading


%number of times structure repreats itself in x, y, and z direction
cell_len_x = 10; cell_len_y = 10; cell_len_z = 10;
period = [len/cell_len_x,bre/cell_len_y,wid/cell_len_z]; 

%name of imported file 
imported_mesh = 'No'; inp = 'import.stl';

%relative density and approximate thickness for uniform density
rd_in = 0.5; app_t = 0.5;

%determines if sheet tpms is evaluated using values of relative density or thickness.
%Yes for relative density, No for thickness
relDenCheck = 'Yes';

%name of tpms surface
lattice_type = 'Primitive Schwartz';    %name of tpms function
custom_function = [];                   % expression of user defined function 

%type of tpms surface
structure_type = 'Sheet';

%determines if the end-caps enclose data values above or below the isovalue
volumeFill = 'below';

%number of ticks in the x, y, and z directions
gridpoints = 150; %numbers should be enclosed in square brackets and seperated by single inverted commas

if strcmp(imported_mesh,'No')       %checks if user want to import a mesh
    %generate gridpoints for x, y, and z
    x = linspace(0, period(1), gridpoints);
    y = linspace(0, period(2), gridpoints);
    z = linspace(0, period(3), gridpoints);
    out_p = [];
elseif strcmp(imported_mesh,'Yes')
    [out_p,x,y,z] = stlinput(inp,gridpoints);
    x = x/cell_len_x; y = y/cell_len_y; z = z/cell_len_z;
else
end

%create meshgrid for x, y, and z
[x1,y1,z1] = meshgrid(x,y,z);

if strcmp(grade,'Yes')  %check if grading is selected
    [isovalue_new,app_t_new] = grading(rd1,rd2,app_t1,app_t2,imported_mesh,inp,len,bre,wid,cell_len_x,cell_len_y,cell_len_z,...
        lattice_type,structure_type, out_p, x1,y1,z1,volumeFill,equation_type,m1,m2,m3,rr,relDenCheck,custom_function,positionX,positionY,positionZ);
    [F,V] = imp_model_func(out_p, x1,y1,z1,cell_len_x,cell_len_y,...
        cell_len_z,app_t_new,isovalue_new,lattice_type,structure_type,volumeFill,custom_function);
else                    %else compute for uniform density
    [isovalue_new,app_t_new] = rel_den(rd_in,imported_mesh,inp,len,bre,wid,cell_len_x,cell_len_y,...
        cell_len_z,app_t,lattice_type,structure_type,out_p,x1,y1,z1,volumeFill,relDenCheck,custom_function);
    [F,V] = imp_model_func(out_p, x1,y1,z1,cell_len_x,cell_len_y,...
        cell_len_z,app_t_new,isovalue_new,lattice_type,structure_type,volumeFill,custom_function);
end

stlwrite2('newfile.stl', F, V); %creates the object on an stl file


