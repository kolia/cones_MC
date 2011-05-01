function X = move_cone( X , x , y , d , PROB , T )

c = X.state(x,y) ;
if ~c  ,  error('trying to move_cone nonexistent cone') ; end

% update contacts and shift_dLLs recursively
action  = @(z,xy)action_LL_shift(z,xy,d,PROB,T) ;
X       = propagate_action(X,x+(y-1)*X.M0,d,action) ;

end