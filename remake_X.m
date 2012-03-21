function Y = remake_X( cone_map, X )

if ~isfield(X,'WW') || ~isfield(X,'contact')
    Y = initialize_X( cone_map.M0, cone_map.M1, ...
                      cone_map.N_colors, cone_map.SS, ...
                      cone_map.cone_params.replusion_radii, ...
                      cone_map.naive_LL, 1, 1) ;

    [x,y,c] = find(X.state) ;
    for i=1:numel(x)
        Y = flip_LL( Y , [x(i) y(i) c(i)] , cone_map , [1 1] ) ;
    end
else
    Y = X ;
end

end