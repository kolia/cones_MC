function X = update_X(X,trial,LL)

type = trial.move{1} ;
x    = trial.move{2} ;
y    = trial.move{3} ;

switch type
    
    case 'delete'
        fprintf('\ndelete %3d %3d',x,y)
        
        X = update_cone_deletion(X,x,y) ;
        X = transitive_closures( X ) ;
        
    case 'add'
        
        % assign new id to added cone
        id              = find( ~X.taken_ids , 1 ) ;
        X.id(x,y)       = id ; 
        X.taken_ids(id) = true ;
        
        % update contacts locally and shift_dLLs recursively
        X               = make_contacts( X , x , y , id , LL ) ;
        X               = transitive_closures( X ) ;
        
        c               = trial.move{4} ;
        X.state(x,y)    = c ;
        X.N_cones       = X.N_cones + 1 ;
        
        fprintf('\nadd    %3d %3d %d',x,y,c)
        
    case 'shift'

        id  = X.id(x,y) ;
        d   = trial.move{4} ;
                
        % update contacts and shift_dLLs recursively
        [X,done]   = propagate_action(X,id,d, @(x,ID)action_shift(x,ID,d,LL)) ;        
        
        check_X( X , LL )
        
        X   = transitive_closures( X ) ;
    
        fprintf('\nshift  %3d %3d %d : %d cones',x,y,d,sum(done))
        
    case 'change color'
        
        c = trial.move{4} ;
        X.state(x,y) = c ;


        fprintf('\ncolor  %3d %3d %d',x,y,c)

end

X.ll = trial.ll ;

end