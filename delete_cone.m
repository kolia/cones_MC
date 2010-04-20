function trial = delete_cone( X , x , y , LL )

old_color = X.state(x,y) ;

if old_color
    trial.ll    = X.ll - LL(x,y,old_color) + X.N_cones_factor ;
    trial.move = {'delete' [x y]} ;    
else
    error('delete_cone cannot delete nonexistent cone...')
end

end