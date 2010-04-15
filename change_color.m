function X = change_color( X , x , y , c , LL )
    
old_color = X.state(x,y) ;

if old_color
    X.ll    = X.ll + LL(x,y,c) - LL(x,y,old_color) ;
    X.state(x,y) = c ;
    
else
    error('change_color tried to change color of nonexisting cone...')
end

end