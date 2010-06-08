function X = move_cone( X , x , y , d , cell_consts , STA_W )

c = X.state(x,y) ;
if ~c  ,  error('trying to move_cone nonexistent cone') ; end

% update contacts and shift_dLLs recursively
action  = @(z,xy)action_LL_shift(z,xy,d,cell_consts,STA_W) ;
X       = propagate_action(X,x+(y-1)*X.M0,d,action) ;

end