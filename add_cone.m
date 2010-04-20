function trial = add_cone( X , x , y , c , LL )

[mask,indices] = place_mask( X.M0 , X.M1 , x , y , X.masks.exclusion ) ;
free    = isempty( find( X.state(indices)>0 , 1) ) ;

if ~free
    trial.ll = -Inf ;
else
    trial.ll = X.ll + LL(x,y,c) - X.N_cones_factor ;    
    trial.move = {'add' x y c} ;
    
end

end