function Y = remake_X( cone_map, X )

Y = X ;

if ~isfield(X,'WW') || ~isfield(X,'contact')
    Y = initialize_X( cone_map.M0, cone_map.M1, ...
                      cone_map.N_colors, cone_map.SS, ...
                      cone_map.cone_params.replusion_radii, ...
                      1, 1) ;

    [x,y,c] = find(X.state) ;
    for i=1:numel(x)
        Y = flip_LL( Y , [x(i) y(i) c(i)] , cone_map , [1 1] ) ;
        Y = update_X({Y}) ;
    end
end

if ~isfield(X, 'excluded')
    Y.excluded  = false(Y.M0,X.M1) ;
    [x,y,c] = find(Y.state) ;
    for i=1:numel(x)
        [~,indices] = not_excluded( Y, x(i), y(i) ) ;
        Y.excluded(indices) = true ;
    end    
end

end