function X = add_cone( X , x , y , c , LL , N_cones_factor , masks )

indices = place_mask( X.M0 , X.M1 , x , y , masks.exclusion ) ;

free    = isempty( find( X.state(indices)>0 , 1) ) ;

if ~free
    X.ll = -Inf ;
else
    X.ll = X.ll + LL(x,y,c) - N_cones_factor ;
    X.state(x,y) = c ;
    X.N_cones = X.N_cones + 1 ;
    
    X.to_update = {'add' x y c} ;
    
end

end