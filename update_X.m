function X = update_X(X,trial,LL)

type = trial.move{1} ;
x    = trial.move{2} ;
y    = trial.move{3} ;
N    = X.N_cones ;
switch type
    
    case 'delete'
%         fprintf('\ndelete %3d %3d',x,y)
        
        X = update_cone_deletion(X,x,y) ;
        X = transitive_closures( X ) ;
        
    case 'add'
        
        % assign new id to added cone
        id              = find( ~X.taken_ids , 1 ) ;
        X.id(x,y)       = id ; 
        X.taken_ids(id) = true ;
        
        % update state color and N_cones
        c               = trial.move{4} ;
        X.state(x,y)    = c ;
        X.N_cones       = X.N_cones + 1 ;

        % update contacts and shifts locally ; all transitive closure
        X               = make_contacts( X , x , y , id , LL ) ;
        X               = transitive_closures( X ) ;
                
%         fprintf('\nadd    %3d %3d %d',x,y,c)
        
    case 'shift'

        id  = X.id(x,y) ;
        d   = trial.move{4} ;
        
        % update contacts and shift_dLLs recursively
        [X,done]   = propagate_action(X,id,d, @(x,ID)action_shift(x,ID,d,LL)) ;
        
%         fprintf('\nshift  %3d %3d %d : %d cones',x,y,d,sum(done))
        
        X   = transitive_closures( X ) ;
       
        
    case 'change color'
        
        c = trial.move{4} ;
        X.state(x,y) = c ;

        for d=1:4
            X.shift_dLL{d}(X.id(x,y)) = local_dLL(LL,x,y,d,X) ;
        end
%         fprintf('\ncolor  %3d %3d %d',x,y,c)

end

X.ll = trial.ll ;
if X.N_cones < N , check_X(LL,X) , end

end