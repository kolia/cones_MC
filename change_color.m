function trial = change_color( X , x , y , c , cell_consts , STA_W )
    
old_color = X.state(x,y) ;

if old_color
    X          = flip_LL( X , [x y c] , cell_consts , STA_W ) ;
    trial.move = {'change color' X x y c} ;
else
    error('change_color tried to change color of nonexisting cone...')
end

end