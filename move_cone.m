function trial = move_cone( X , x , y , d , LL )

c = X.state(x,y) ;
if ~c  ,  error('trying to move_cone nonexistent cone') ; end

move    = [x y] + X.masks.shift{d} ;

% moving cone outside of borders: just delete cone
if  ~( move(1)  &&  move(2)  &&  X.M0>=move(1)  &&  X.M1>=move(2) )
    trial = delete_cone( X , x , y , LL ) ;
    
% regular move of cone x,y
else
    trial.ll = X.ll + X.shift_dLL{d}( X.id(x,y) ) ;
    trial.move = {'shift' x y d} ;
end



end