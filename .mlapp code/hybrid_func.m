function [F,V] = hybrid_func(app,rd1,app_t1,imported_mesh,inp,len,bre,wid,out_p,x1,y1,z1,cell_len_x,cell_len_y, cell_len_z,lattice_type,lattice_type2,structure_type,volumeFill, ...
    custom_function,custom_function2,k,hybridEquation,a,b,c,rr,relDenCheck,positionX,positionY,positionZ)
 %check if app forced termination
drawnow;
if app.stopFlag
    return
end

switch volumeFill
    case "Fill Above Isovalue"
        volumeFill = 'above';
    case "Fill Below Isovalue"
        volumeFill = 'below';
    otherwise 
end

f1 = tpms_function(app,x1,y1,z1,lattice_type,custom_function);                  %gets the value of f for corresponding tpms structure
f2 = tpms_function(app,x1,y1,z1,lattice_type2,custom_function2);
f1(out_p == 0) = 1e10; f2(out_p == 0) = 1e10;                               %redefine level set field based in case of arbitrary geometry

x_c = positionX; y_c = positionY; z_c = positionZ;
xx1 = x1-x_c/cell_len_x; yy1 = y1-y_c/cell_len_y; zz1 = z1-z_c/cell_len_z;     %define the reference frame
xyz_max = max([max(xx1(:)) max(yy1(:)) max(zz1(:))]);                          %obtain the maximum dimension
xx1 = xx1/xyz_max; yy1 = yy1/xyz_max; zz1 = zz1/xyz_max;

if strcmp(hybridEquation,'Linear')                                           %for linear
    g = a*(xx1) +b*(yy1) + c*(zz1);
elseif strcmp(hybridEquation,'Quadratic')                                    %quadratic
    g = a*(xx1.^2) +b*(yy1.^2) + c*(zz1.^2);
elseif strcmp(hybridEquation,'Cubic')                                        %cubic
    g = a*(xx1.^3) +b*(yy1.^3) + c*(zz1.^3);
elseif strcmp(hybridEquation,'Circular')                                     %for circular
    rr = rr/(max([cell_len_x cell_len_y cell_len_z]).*xyz_max);
    g = sqrt(a.*xx1.^2+b.*yy1.^2+c.*zz1.^2)-rr;
elseif strcmp(hybridEquation,'Exponential')
    g = exp(a*xx1) + exp(b*yy1) + exp(c*zz1);
elseif strcmp(hybridEquation,'Trigonometric - Sin')
    g = sin(a*xx1) + sin(b*yy1) + sin(c*zz1);
elseif strcmp(hybridEquation,'Trigonometric - Cosine')
    g = cos(a*xx1) + cos(b*yy1) + cos(c*zz1);
end

switch structure_type
    case "Sheet"

        if relDenCheck == true
            [~,at1] = rel_den(app,rd1,imported_mesh,inp,len,bre,wid,cell_len_x,cell_len_y,cell_len_z,... %uses relative density to determine thickness of grading
                [],lattice_type,structure_type, out_p, x1,y1,z1,volumeFill,relDenCheck,custom_function);
        else
            at1 = app_t1;                                                             %else uses thickness value
        end

        isovalue = at1/2;

        %hybridization formula
        s = 1/(1 + exp(k.*g));
        p = s.*f1 + (1-s).*f2;
        f = p;

        h = (f-isovalue).*-(f+isovalue);

        drawnow;
        if app.stopFlag
            return
        end

        [face1, vertex1] = isosurface(x1,y1,z1,h,0);    %creates isosurface

        drawnow;
        if app.stopFlag
            return
        end

        [face2, vertex2] = isocaps(x1,y1,z1,h,0);       %creates endcaps

    case "Solid"

        [isovalue1,~] = rel_den(app,rd1,imported_mesh,inp,len,bre,wid,cell_len_x,cell_len_y,cell_len_z,... %uses relative density to determine thickness of grading
            [],lattice_type,structure_type, out_p, x1,y1,z1,volumeFill,relDenCheck,custom_function);
        [isovalue2,~] = rel_den(app,rd1,imported_mesh,inp,len,bre,wid,cell_len_x,cell_len_y,cell_len_z,...
            [],lattice_type2,structure_type, out_p, x1,y1,z1,volumeFill,relDenCheck,custom_function);

        h1 = f1+isovalue1;
        h2 = f2+isovalue2;

        s = 1/(1 + exp(k.*g));
        p = s.*h1 + (1-s).*h2;
        h = p;

        drawnow;
        if app.stopFlag
            return
        end

        [face1, vertex1] = isosurface(x1,y1,z1,h,0);

        drawnow;
        if app.stopFlag
            return
        end

        [face2, vertex2] = isocaps(x1,y1,z1,h,0, volumeFill );
end

drawnow;
if app.stopFlag
    return
end

F = [face1; length(vertex1(:,1)) + face2]; %concatenates the face arrays which represent the faces of the resulting Solid
V = [vertex1; vertex2]; %concatenates the vertex arrays which represent the vertices of the resulting Solid
V = [V(:,1)*cell_len_x, V(:,2)*cell_len_y, V(:,3)*cell_len_z];
