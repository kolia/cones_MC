function Y = remake_X( cone_map, X )

if ~isfield(X,'WW')
    Y = cone_map.initX ;

    [x,y,c] = find(X.state) ;
    for i=1:numel(x)
        Y = flip_LL( Y , [x(i) y(i) c(i)] , cone_map , [1 1] ) ;
    end
else
    Y = X ;
end

end