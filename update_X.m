function X = update_X(X,trial,LL)

type = trial.move{1} ;
Y    = trial.move{2} ;
x    = trial.move{3} ;
y    = trial.move{4} ;

X.overlaps  = Y.overlaps ;
X.invWW     = Y.invWW ;       % X.WW = Y.WW ;
X.positions = Y.positions ;
X.colors    = Y.colors ;
X.state     = Y.state ;
X.ll        = Y.ll ;

switch type
    
    case 'delete'
%         fprintf('\ndelete %3d %3d',x,y)
        X   = update_cone_deletion(X,x,y) ;
        
    case 'add'
                                                        % WHERE DO I UPDATE X.state ???
        if X.N_cones >= X.maxcones                      % WHICH CONE DO I DELETE    ??? !!!
            [minLL,mid] = min(X.localLL(X.localLL>0)) ;
            [mx   ,my ] = find(X.id==mid,1) ;
            X = update_cone_deletion(X,mx,my) ;
        end
        
        c   = trial.move{5} ;
        X   = update_cone_addition(X,x,y,c) ;
                
%         fprintf('\nadd    %3d %3d %d',x,y,c)
        
    case 'shift'

        d   = trial.move{5} ;
        
        % update contacts and shift_dLLs recursively
        X   = propagate_action(X,x,y,d,@(X,x,y)action_shift(X,x,y,d)) ;
        
%         fprintf('\nshift  %3d %3d %d : %d cones',x,y,d,sum(done))
        
               
    case 'change color'
        
%         fprintf('\ncolor  %3d %3d %d',x,y,c)

end

X   = transitive_closures( X ) ;

check_X(LL,X)

end