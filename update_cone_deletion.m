function X = update_cone_deletion(X,x,y)

% free cone's id
id              = X.id(x,y) ;
X.id(x,y)       = 0 ;
X.taken_ids(id) = false ;
X.state(x,y)    = 0 ;
X.N_cones       = X.N_cones - 1 ;

% update contacts locally and shift_dLLs recursively
for d=1:4
    X           = delete_contact( X , id , d ) ;
end


end