function X = update_shift(X,x,y,d)

% check that move is possible without moving other cones
[mask,indices]    = place_mask( X.M0 , X.M1 , x , y , X.masks.cardinal{d} ) ;
contact = find( X.state(indices) ) ;
free    = isempty( contact ) ;

force   = 0 ;
if ~free
    for i = contact
        X = update_shift(X,mask(i,1),mask(i,2),d) ;
        force = force +  ;
    end    
end

move    = [x y] + X.masks.shift{d} ;
X.state(move(1),move(2)) = X.state(x,y) ;
X.state(x,y) = 0 ;

end