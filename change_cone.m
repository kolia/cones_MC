function trial = change_cone( X , a , cell_consts , STA_W )

for i=1:size(a,1)
    x = a(i,1) ;
    y = a(i,2) ;
    c = a(i,3) ;
    free = 1 ;
    
    if c && ~X.state(x,y)        
        [mask,indices] = place_mask( X.M0 , X.M1 , x , y , X.masks.exclusion ) ;
        free    = isempty( find( X.state(indices)>0 , 1) ) ;
    end

    if free
        X    = flip_LL( X , [x y c] , cell_consts , STA_W ) ;
    else
        X.ll = -Inf ;
    end
    
end

trial.move = {'generate' X a} ;

end