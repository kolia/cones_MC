function X = action_LL_shift(X,ox,oy,d,cell_consts,STA_W)

% new cone location
x   = ox+X.masks.shift{d}(1) ;
y   = oy+X.masks.shift{d}(2) ;
c   = X.state(ox,oy) ;    
id  = X.id(ox,oy) ;

% fprintf('\n\nMOVING cone %d at %d,%d',id,ox,oy)

% delete old cone
X = flip_LL( X , [ox oy 0] , cell_consts , STA_W ) ;

X.id(ox,oy)     = 0  ;
X.state(ox,oy)  = 0 ;
% delete old contacts
for dd=1:4
    X = delete_contact(X,id,dd) ;
end

% if new cone location is within bounds, add new cone
if x>0 && x<=X.M0 && y>0 && y<=X.M1
    X = flip_LL( X , [x y c] , cell_consts , STA_W ) ;
        
    % shift cone state and id
    X.state(x,y)    = c ;
    X.id( x, y)     = id ;
    
    % make new contacts
    X = make_contacts(X,x,y,id) ;
    
else
    X.taken_ids(id) = false ;
    
end

end