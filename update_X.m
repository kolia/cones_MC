function X = update_X(X,trial,LL)

type = trial.move{1} ;
x    = trial.move{2} ;
y    = trial.move{3} ;

X.ll = trial.ll ;

switch type
    
    case 'delete'
%         fprintf('\ndelete %3d %3d',x,y)
        
        X   = update_cone_deletion(X,x,y) ;
        X   = transitive_closures( X ) ;
        
    case 'add'
        
        if X.N_cones >= X.maxcones
            [minLL,mid] = min(X.localLL(X.localLL>0)) ;
            [mx   ,my ] = find(X.id==mid,1) ;
            X = update_cone_deletion(X,mx,my) ;
        end
        
        c   = trial.move{4} ;
        X   = update_cone_addition(X,x,y,c,LL) ;        
        
        % all transitive closures
        X   = transitive_closures( X ) ;
                
%         fprintf('\nadd    %3d %3d %d',x,y,c)
        
    case 'shift'

        id  = X.id(x,y) ;
        d   = trial.move{4} ;
        
        % update contacts and shift_dLLs recursively
        X   = propagate_action(X,id,d, @(x,ID)action_shift(x,ID,d,LL)) ;
        
%         fprintf('\nshift  %3d %3d %d : %d cones',x,y,d,sum(done))
        
        X   = transitive_closures( X ) ;
       
        
    case 'change color'
        
        c   = trial.move{4} ;
        X.state(x,y) = c ;
        id  = X.id(x,y) ;

        for d=1:4
            X.shift_dLL{d}(id) = local_dLL(LL,x,y,d,X) ;
        end
        X.localLL(id)   = LL(x,y,c) ;
%         fprintf('\ncolor  %3d %3d %d',x,y,c)

end

% if X.N_cones < N , check_X(LL,X) , end
check_X(LL,X)

end