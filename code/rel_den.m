function [isovalue_new,app_t_new] = rel_den(rd_in,imported_mesh,inp,len,bre,wid,cell_len_x,cell_len_y,...
    cell_len_z,app_t,lattice_type,structure_type,out_p,x1,y1,z1,volumeFill,relDenCheck,custom_function)

if strcmp(structure_type,'Solid')                                                           % Check if Solid
    isovalue = [-1 1]; rd = zeros(size(isovalue));                                          % Initialize isovalue for two extreme values -1 and 1
    for i = 1:length(isovalue)                                                              % for-loop for the two isovalues chosen
        [F,V] = imp_model_func(out_p,x1,y1,z1,cell_len_x,cell_len_y,...
            cell_len_z,app_t,isovalue(i),lattice_type,structure_type,volumeFill,custom_function);
        if strcmp(imported_mesh,'Yes')                                                      % Check if the user is using an imported mesh domain
            nd = stlread(inp);                                                              % read the imported mesh domain
            ndv = nd.Points; ndf = nd.ConnectivityList;                                     % get face and vertices of imported mesh domain
            v_ini = stlVolume(ndv',ndf');                                                   % compute the volume of imported mesh domain
        else
            v_ini = len*bre*wid;                                                            % compute volume of default domain
        end
        rd(i) = abs(stlVolume(V',F'))/v_ini;                                                % calculate the relative density of TPMS based on a particular isovalue
    end
    isovalue_new = (isovalue(2)-isovalue(1))*((rd_in-rd(1))/(rd(2)-rd(1)))+isovalue(1);     % estimate the new isovalue by interpolating between the two initial extreme isovalues
    app_t_new = [];
else
    if strcmp(relDenCheck, "Yes")                                                           %check if new approximate thickness is calculated using relative density
        app_t= [0.2 2];
        rd = zeros(size(app_t));
        for i = 1:length(app_t)
            [F,V] = imp_model_func(out_p,x1,y1,z1,cell_len_x,cell_len_y,...
                cell_len_z,app_t(i),[],lattice_type,structure_type,volumeFill,custom_function);
            if strcmp(imported_mesh,'Yes')
                nd = stlread(inp);
                ndv = nd.Points; ndf = nd.ConnectivityList;
                v_ini = stlVolume(ndv',ndf');
            else
                v_ini = len*bre*wid;
            end
            rd(i) = abs(stlVolume(V',F'))/v_ini;
        end
        app_t_new = (app_t(2)-app_t(1))*((rd_in-rd(1))/(rd(2)-rd(1)))+app_t(1);
        isovalue_new = [];
    else
        app_t_new = app_t;
        isovalue_new = [];
    end
end

