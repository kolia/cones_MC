function X = update_cone_addition(X,x,y,c,LL)

% assign new id to added cone
id              = find( ~X.taken_ids , 1 ) ;

X.id(x,y)       = id ;
X.taken_ids(id) = true ;

% update state color and N_cones
X.state(x,y)    = c ;
X.N_cones       = X.N_cones + 1 ;

% update contacts and shifts locally
X               = make_contacts( X , x , y , id , LL ) ;

end