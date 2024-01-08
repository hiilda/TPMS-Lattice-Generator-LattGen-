function [isovalue_new,app_t_new] = grading(app,rd1,rd2,app_t1,app_t2,imported_mesh,inp,len,bre,wid,cell_len_x,cell_len_y,cell_len_z,...
    lattice_type,structure_type, out_p, x1,y1,z1,volumeFill,equation_type,m1,m2,m3,rr,relDenCheck,custom_function,positionX,positionY,positionZ)
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

x_c = positionX; y_c = positionY; z_c = positionZ;      %reference frame

switch structure_type
    case "Sheet"
        if relDenCheck == true
            [~,at1] = rel_den(app,rd1,imported_mesh,inp,len,bre,wid,cell_len_x,cell_len_y,cell_len_z,... %uses relative density to determine thickness of grading
                [],lattice_type,structure_type, out_p, x1,y1,z1,volumeFill,relDenCheck,custom_function);
            [~,at2] = rel_den(app,rd2,imported_mesh,inp,len,bre,wid,cell_len_x,cell_len_y,cell_len_z,...
                [],lattice_type,structure_type, out_p, x1,y1,z1,volumeFill,relDenCheck,custom_function);
        else
            at1 = app_t1; at2 = app_t2;                                                             %else uses thickness value 
        end

        xx = x1-x_c/cell_len_x; yy = y1-y_c/cell_len_y; zz = z1-z_c/cell_len_z;     %define the reference frame
        xyz_max = max([max(xx(:)) max(yy(:)) max(zz(:))]);                          %obtain the maximum dimension
        xx = xx/xyz_max; yy = yy/xyz_max; zz = zz/xyz_max;  %xx1,yy1,zz1            %normalize with the maximum dimension

        if strcmp(equation_type,'Linear')                                           %for linear grading
            is_f = abs(m1*xx+m2*yy+m3*zz);
        elseif strcmp(equation_type,'Quadratic')                                    %quadratic grading
            is_f = abs(m1*xx.^2+m2*yy.^2+m3*zz.^2);
        elseif strcmp(equation_type,'Cubic')                                        %cubic grading
            is_f = abs(m1*xx.^3+m2*yy.^3+m3*zz.^3);
        elseif strcmp(equation_type,'Circular')                                     %for circular grading patterns
            rr = rr/max(max([cell_len_x cell_len_y cell_len_z]).*xyz_max);
            is_f = abs(sqrt(m1.*xx.^2+m2.*yy.^2+m3.*zz.^2)-rr);
        elseif strcmp(equation_type,'Exponential')
            is_f = abs(exp(m1*xx)+exp(m2*yy)+exp(m3*zz));
        elseif strcmp(equation_type, 'Trigonometric - Sin')
            is_f = abs(sin(m1*xx)+sin(m2*yy)+sin(m3*zz));
        elseif strcmp(equation_type, 'Trigonometric - Cosine')
            is_f = abs(cos(m1*xx)+cos(m2*yy)+cos(m3*zz));
        end
        app_t_new = (at1+is_f*(at2-at1));                                           %interpolation function to obtain new and unique isovalues
        isovalue_new = (at1+is_f*(at2-at1))/2;

    case "Solid"

        [is1,~] = rel_den(app,rd1,imported_mesh,inp,len,bre,wid,cell_len_x,cell_len_y,cell_len_z,...
            [],lattice_type,structure_type, out_p, x1,y1,z1,volumeFill,relDenCheck,custom_function);
        [is2,~] = rel_den(app,rd2,imported_mesh,inp,len,bre,wid,cell_len_x,cell_len_y,cell_len_z,...
            [],lattice_type,structure_type, out_p, x1,y1,z1,volumeFill,relDenCheck,custom_function);

        xx = x1-x_c/cell_len_x; yy = y1-y_c/cell_len_y; zz = z1-z_c/cell_len_z;
        xyz_max = max([max(xx(:)) max(yy(:)) max(zz(:))]);
        xx = xx/xyz_max; yy = yy/xyz_max; zz = zz/xyz_max;
        if strcmp(equation_type,'Linear')
            is_f = abs(m1*xx+m2*yy+m3*zz);
        elseif strcmp(equation_type,'Quadratic')
            is_f = abs(m1*xx.^2+m2*yy.^2+m3*zz.^2);
        elseif strcmp(equation_type,'Cubic')
            is_f = abs(m1*xx.^3+m2*yy.^3+m3*zz.^3);
        elseif strcmp(equation_type,'Circular') 
            rr = rr/max([len bre wid]);
            is_f = abs(sqrt(m1.*xx.^2+m2.*yy.^2+m3.*zz.^2)-rr);
        elseif strcmp(equation_type,'Exponential')
            is_f = abs(exp(m1*xx)+exp(m2*yy)+exp(m3*zz));
        elseif strcmp(equation_type, 'Trigonometric - Sin')
            is_f = abs(sin(m1*xx)+sin(m2*yy)+sin(m3*zz));
        elseif strcmp(equation_type, 'Trigonometric - Cosine')
            is_f = abs(cos(m1*xx)+cos(m2*yy)+cos(m3*zz));
        end
        isovalue_new = is1+is_f*(is2-is1);
        app_t_new = [];
end
end
