function trial = move_cone( X , x , y , d , cell_consts , STA_W )

c = X.state(x,y) ;
if ~c  ,  error('trying to move_cone nonexistent cone') ; end

move    = [x y] + X.masks.shift{d} ;

% moving cone outside of borders: just delete cone
if  ~( move(1)  &&  move(2)  &&  X.M0>=move(1)  &&  X.M1>=move(2) )
    trial = delete_cone( X , x , y , cell_consts , STA_W ) ;
    
% regular move of cone x,y
else
    
    % update contacts and shift_dLLs recursively
    action  = @(z,x,y)action_LL_shift(z,x,y,d,cell_consts,STA_W) ;
    X       = propagate_action(X,x,y,d,action) ;

    trial.move = {'shift' X [x y d]} ;
end

end