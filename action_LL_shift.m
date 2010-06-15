function X = action_LL_shift(X,oxy,d,PROB)

ox  = 1 + mod(oxy-1,X.M0) ;
oy  = 1 + floor((oxy-1)/X.M0) ;

% new cone location
x   = ox+X.masks.shift{d}(1) ;
y   = oy+X.masks.shift{d}(2) ;
c   = X.state(ox,oy) ;    

% fprintf('\n\nMOVING cone %d at %d,%d',id,ox,oy)

% delete old cone
X = flip_LL( X , [ox oy 0] , PROB ) ;

% if new cone location is within bounds, add new cone
if x>0 && x<=X.M0 && y>0 && y<=X.M1
    X = flip_LL( X , [x y c] , PROB ) ;
end

end