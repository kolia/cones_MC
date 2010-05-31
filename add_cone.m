function trial = add_cone( X , x , y , c , cell_consts , STA_W )

[mask,indices] = place_mask( X.M0 , X.M1 , x , y , X.masks.exclusion ) ;
free    = isempty( find( X.state(indices)>0 , 1) ) ;

if ~free
    X.ll = -Inf ;
else
    
    % delete least probable cone if maximum number of cones reached
    if X.N_cones == X.maxcones
        [xmin,ymin] = find( X.id == find(X.localLL==min(X.localLL)) ) ;
        X = flip_LL( X , [xmin ymin 0] , cell_consts , STA_W ) ;
    end
    
    X     = flip_LL( X , [x y c] , cell_consts , STA_W ) ;
    
    trial.move = {'add' X x y c} ;
    
end

end