function X = action_LL_shift(X,ox,oy,d,cell_consts,STA_W)

% new cone location
x = ox+X.masks.shift{d}(1) ;
y = oy+X.masks.shift{d}(2) ;
c = X.state(ox,oy) ;

% delete old cone
X = flip_LL( X , [ox oy 0] , cell_consts , STA_W ) ;
    
% if new cone location is within bounds, add new cone
if x>0 && x<=X.M0 && y>0 && y<=X.M1
    X = flip_LL( X , [x y c] , cell_consts , STA_W ) ;
end

end