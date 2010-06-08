function X = update_X(X,trial)

check_X(X)

type = trial.move{1} ;
Y    = trial.move{2} ;
a    = trial.move{3} ;

X.overlaps  = Y.overlaps ;
X.invWW     = Y.invWW ;       % X.WW = Y.WW ;
X.positions = Y.positions ;
X.colors    = Y.colors ;
X.state     = Y.state ;
X.ll        = Y.ll ;
X.N_cones   = Y.N_cones ;

n_changes = size(a,1) ;

switch type
    
    case 'generate'
        for i=1:n_changes
            x = a(i,1) ;
            y = a(i,2) ;
            c = a(i,3) ;

            % add cone
            if c
                if n_changes == 1  % unless color change, assign new id
                    id              = find( ~X.taken_ids , 1 ) ;
                    X.id(x,y)       = id ;
                    X.taken_ids(id) = true ;
                    X               = make_contacts( X , x , y , id ) ;
                end

            % delete cone
            else
                if n_changes == 1  % unless color change, free id
                    id              = X.id(x,y) ;
                    X.id(x,y)       = 0 ;
                    X.taken_ids(id) = false ;
                    
                    % update contacts locally
                    for d=1:4
                        X           = delete_contact( X , id , d ) ;
                    end
                end            
            end                
        end        
        fprintf('\nchange %3d %3d %d     ll: %f',x,y,c,X.ll)
        
    case 'shift'% update contacts recursively
        d = a(1,3) ;
        id= X.id(a(1,1),a(1,2)) ;
%         X = propagate_action(X,a(1,1),a(1,2),id,d,@(X,x,y)action_shift(X,x,y,d)) ;
        

        x = a(1,1) ;
        y = a(1,2) ;
        fprintf('\nshift  %3d %3d %d',x,y,d)

end

X = transitive_closures( X ) ;
check_X(X)

end