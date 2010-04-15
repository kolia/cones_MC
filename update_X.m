function X = update_X(X)

if isfield(X.to_update)

    type = X.to_update{1} ;
    x    = X.to_update{2} ;
    y    = X.to_update{3} ;

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

        case 'add'

            id          = find(X.taken_ids , 1) ;
            X.id(x,y)   = id ;

            X.taken_ids(id) = true ;

            for d=1:4
                inds =  x + X.masks.cardinal{d}(:,1)      + ...
                       (y + X.masks.cardinal{d}(:,2) - 1) * X.M0 ;
                X.contact{d}( inds ) = X.state( inds ) ;
            end

            
        case 'move'

    end

    X = rmfield(X , 'to_update') ;
end

end