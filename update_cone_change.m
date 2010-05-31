function X = update_cone_change(X,x,y,c)


if c% add cone

    % assign new id to added cone
    id              = find( ~X.taken_ids , 1 ) ;
    
    X.id(x,y)       = id ;
    X.taken_ids(id) = true ;
    
    % update N_cones
    X.N_cones       = X.N_cones + 1 ;
    
    % update contacts and shifts locally
    X               = make_contacts( X , x , y , id ) ;
    

else% delete cone
    
    % free cone's id
    id              = X.id(x,y) ;
    X.id(x,y)       = 0 ;
    X.taken_ids(id) = false ;
    X.N_cones       = X.N_cones - 1 ;
    
    % update contacts locally and shift_dLLs recursively
    for d=1:4
        X           = delete_contact( X , id , d ) ;
    end
    
end
    
end