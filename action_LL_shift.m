function X = action_LL_shift(X,oxy,d,PROB,T)

ox  = 1 + mod(oxy-1,X.M0) ;
oy  = 1 + floor((oxy-1)/X.M0) ;

% new cone location
x   = ox+X.masks{1,1}.shift{d}(1) ;
y   = oy+X.masks{1,1}.shift{d}(2) ;
c   = X.state(ox,oy) ;    

% fprintf('\n\nMOVING cone %d at %d,%d',id,ox,oy)

% delete old cone
X = flip_LL( X , [ox oy 0] , PROB , T ) ;
% check_X(X)

% if new cone location is within bounds, add new cone
if x>0 && x<=X.M0 && y>0 && y<=X.M1
    X = flip_LL( X , [x y c] , PROB , T ) ;
%     check_X(X)
end

end