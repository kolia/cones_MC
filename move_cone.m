function trial = move_cone( X , x , y , d , cell_consts , STA_W )

oll = X.ll ;

c = X.state(x,y) ;
if ~c  ,  error('trying to move_cone nonexistent cone') ; end

% update contacts and shift_dLLs recursively
action  = @(z,x,y)action_LL_shift(z,x,y,d,cell_consts,STA_W) ;
X       = propagate_action(X,x,y,X.id(x,y),d,action) ;

% fprintf('\nmove dll %f',X.ll-oll)

trial.move = {'shift' X [x y d]} ;
trial.ll   = X.ll ;

end