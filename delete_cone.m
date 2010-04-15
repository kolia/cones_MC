function X = delete_cone( X , x , y , LL , N_cones_factor )

old_color = X.state(x,y) ;

if old_color
    X.ll        = X.ll - LL(x,y,old_color) + N_cones_factor ;
    X.state(x,y)= 0 ;
    X.N_cones   = X.N_cones - 1 ;

    X.to_update = {'delete' x y 0} ;
    
else
    error('delete_cone cannot delete nonexistent cone...')
end

end