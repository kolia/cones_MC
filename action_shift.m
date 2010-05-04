function X = action_shift(X,id,d,LL)

% old cone location
[ox,oy] = find( X.id == id ) ;
if length(ox) ~= 1
    error('X.id should contain exactly 1 entry with this id.') ;
end

% new cone location
x = ox+X.masks.shift{d}(1) ;
y = oy+X.masks.shift{d}(2) ;


% Y = X ;

% if new cone location is within bounds
if x>0 && x<=X.M0 && y>0 && y<=X.M1

    % delete old contacts
    for dd=1:4
        X = delete_contact(X,id,dd) ;
    end
    
    % shift cone state and id
    X.id( x, y) = id ;
    X.id(ox,oy) = 0  ;
    X.state( x, y) = X.state(ox,oy) ;
    X.state(ox,oy) = 0 ;
    
    % make new contacts and get local dLLs
    X = make_contacts(X,x,y,id,LL) ;
    
% if new cone location is out of bounds
else
    X = update_cone_deletion(X,ox,oy) ;
end


end