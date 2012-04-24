function cone_map = remake_initX( cone_map, hot )
% remake cone_map.initX if it was deleted

if ~isfield(cone_map,'initX')
    if nargin<2, hot = false ; end

    %% initial empty X
    cone_map.initX = initialize_X( cone_map.M0, cone_map.M1, ...
                                   cone_map.N_colors, cone_map.SS, ...
                                   cone_map.cone_params.replusion_radii, ...
                                   1, 1) ;

    if hot
        % rerun 'hot' greedy algorithm, to initialize X replay
        greed_hot = GREEDY(cone_map, 'hot') ;
        cone_map.initX = greed_hot.X ;
    end

    %% transfer all info from cone_map to cone_map.initX
    cone_map.initX = transfer_info( cone_map, cone_map.initX ) ;
end

end