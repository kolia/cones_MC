function trial = delete_cone( X , x , y , cell_consts , STA_W )

old_color = X.state(x,y) ;

if old_color
    X          = flip_LL( X , [x y 0] , cell_consts , STA_W ) ;
    trial.move = {'delete' X x y} ;
else
    error('delete_cone cannot delete nonexistent cone...')
end

end