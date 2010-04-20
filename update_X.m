function X = update_X(X,trial)

type = trial.move{1} ;
x    = trial.move{2} ;
y    = trial.move{3} ;

switch type
    
    case 'delete'
        
        id          = X.id(x,y) ;
        X.id(x,y)   = 0 ;
        
        X.taken_ids(id) = false ;
        
        % apply the following pseudocode to all 4 cardinal directions
        % N = X.contact.N(id,:) ;
        % X.contact.S(N,id) = false ;
        for d=1:4
            there = X.contact{d}(id,:) ;
            X.contact{1+mod(d+1,4)}(there,id) = false ;
        end
        
        X.state(x,y)= 0 ;
        X.N_cones   = X.N_cones - 1 ;
        
        
    case 'add'
        
        id          = find( ~X.taken_ids , 1 ) ;
        X.id(x,y)   = id ;
        
        X.taken_ids(id) = true ;
        
        for d=1:4
            mask = place_mask( X.M0 , X.M1 , x , y , X.masks.cardinal{d}) ;
            inds =  mask(:,1) + (mask(:,2) - 1) * X.M0 ;
            X.contact{d}( inds ) = X.state( inds ) ;
        end
        
        c = trial.move{4} ;
        X.state(x,y) = c ;
        X.N_cones = X.N_cones + 1 ;
        
        
    case 'shift'
        
        d = trial.move{4} ;
        X = update_shift(X,x,y,d) ;
        
    case 'change color'
        
        c = trial.move{4} ;
        X.state(x,y) = c ;
end

X.ll = trial.ll ;

end