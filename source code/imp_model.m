clc
clear all
tic

% cuboid domain dimensions
len = 60; bre = 60; wid = 60;

%relative density and approximate thickness for uniform density
rd_in = 0.5; app_t = 0.5;

%hybrid variables 
hybrid = 'No'; hybridEquation = "Linear"; a=1; b=0; c=0; %coefficient for variables in equation
k = 10;  %width of transition

%grading variables
grade = 'Yes'; equation_type = 'Linear'; %user selects type of grading to use 
m1 = 0; m2 = 0; m3 = 1; % coefficients for variables in equations

rr = 2; % radius of circle if needed by user
rd1 = 0.1; rd2 = 0.5; app_t1 = 0.3; app_t2 = 2; % relative densities and approximate thicknesses
positionX = len/2; positionY = bre/2; positionZ = wid; %reference position 

%number of times structure repreats itself in x, y, and z direction
cell_len_x = 10; cell_len_y = 10; cell_len_z = 10;
period = [len/cell_len_x,bre/cell_len_y,wid/cell_len_z]; 

%name of imported file 
imported_mesh = 'No'; inp = '.stl';

%determines if sheet tpms is evaluated using values of relative density or thickness.
%true for relative density, false for thickness
relDenCheck = true;

%name of tpms surface
lattice_type = 'Primitive Schwartz Surface';   %name of tpms function
lattice_type2 = 'Gyroid Surface';               %second function for hybridization
custom_function = [];                   % expression of user defined function 
custom_function2 = [];                  %second function for hybridization

%type of tpms surface
structure_type = 'Sheet';

%determines if the end-caps enclose data values above or below the isovalue
volumeFill = 'above';

%number of ticks in the x, y, and z directions
gridpoints = 100; 

if strcmp(imported_mesh,'No')       %checks if user wants to import a mesh
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
elseif strcmp(hybrid, 'Yes')    %check if hybridization is selected
    [F,V] = hybrid_func(rd1,app_t1,imported_mesh,inp,len,bre,wid,out_p,x1,y1,z1,cell_len_x,cell_len_y, cell_len_z,lattice_type,lattice_type2,structure_type,volumeFill, ...
        custom_function,custom_function2,k,hybridEquation,a,b,c,rr,relDenCheck,positionX,positionY,positionZ);
else                    %else compute for uniform density
    [isovalue_new,app_t_new] = rel_den(rd_in,imported_mesh,inp,len,bre,wid,cell_len_x,cell_len_y,...
        cell_len_z,app_t,lattice_type,structure_type,out_p,x1,y1,z1,volumeFill,relDenCheck,custom_function);
    [F,V] = imp_model_func(out_p, x1,y1,z1,cell_len_x,cell_len_y,...
        cell_len_z,app_t_new,isovalue_new,lattice_type,structure_type,volumeFill,custom_function);
end

stlwrite2('newfile.stl', F, V); %creates the object on an stl file
toc;

