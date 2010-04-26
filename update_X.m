function X = update_X(X,trial)

type = trial.move{1} ;
x    = trial.move{2} ;
y    = trial.move{3} ;

switch type
    
    case 'delete'
        
        % free cone's id
        id          = X.id(x,y) ;
        X.id(x,y)   = 0 ;
        X.taken_ids(id) = false ;
        
        % update contacts locally and shift_dLLs recursively
        for d=1:4
            there = X.contact{d}(id,:) ;
            X.contact{d}(id,there)              = false ;
            X.contact{1+mod(d+1,4)}(there,id)   = false ;
            dLL   = -X.shift_dLL{d}(id) ;
            X = propagate_action(X,id,d, @(x,ID)action_dLL(x,ID,d,dLL)) ;
        end
        
        X.state(x,y)= 0 ;
        X.N_cones   = X.N_cones - 1 ;
        
        
    case 'add'
        
        % assign new id to added cone
        id              = find( ~X.taken_ids , 1 ) ;
        X.id(x,y)       = id ;
        X.taken_ids(id) = true ;
        
        % update contacts locally and shift_dLLs recursively
        for d=1:4
            [mask,inds] = place_mask( X.M0 , X.M1 , x , y , X.masks.cardinal{d}) ;
            X.contact{d}( id , X.id(inds) )             = true ;
            X.contact{1+mod(d+1,4)}( X.id(inds) , id )  = true ;
            
            dLL= local_dLL(X,id,d) ;
            od = 1+mod(d+1,4) ;
            X  = propagate_action(X,id,od, @(x,ID)action_dLL(x,ID,d,dLL)) ;
            
        end

        c = trial.move{4} ;
        X.state(x,y) = c ;
        X.N_cones = X.N_cones + 1 ;
        
        
    case 'shift'

        id  = X.id(x,y) ;
        d   = trial.move{4} ;
        od  = 1+mod(d+1,4) ;
        
        % update contacts locally: sever contacts in opposite direction
        there = X.contact{od}(id,:) ;
        X.contact{od}(id,there)   = false ;
        X.contact{ d}(there,id)   = false ;
        
        % update contacts and shift_dLLs recursively
        X  = propagate_action(X,id,d, @(x,ID)action_shift(x,ID,d)) ;

        X.state(x+X.masks.shift{d}(1),y+X.masks.shift{d}(2)) = X.state(x,y) ;
        X.state(x,y) = 0 ;

    
    case 'change color'
        
        c = trial.move{4} ;
        X.state(x,y) = c ;
end

X.ll = trial.ll ;

end