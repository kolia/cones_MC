function X = update_X(X,trial,LL)

type = trial.move{1} ;
Y    = trial.move{2} ;
a    = trial.move{3} ;

X.overlaps  = Y.overlaps ;
X.invWW     = Y.invWW ;       % X.WW = Y.WW ;
X.positions = Y.positions ;
X.colors    = Y.colors ;
X.state     = Y.state ;
X.ll        = Y.ll ;

switch type
    
    case 'generate'
        for i=1:size(a,1)
            x = a(i,1) ;
            y = a(i,2) ;
            c = a(i,3) ;

            if c && ~X.state(x,y)
                X   = update_cone_addition(X,x,y,c) ;
            elseif ~c
                X   = update_cone_deletion(X,x,y) ;
            end                
        end        
%         fprintf('\nchange %3d %3d %d',x,y,c)
        
    case 'shift'
        d = a(1,3) ;
        
        % update contacts and shift_dLLs recursively
        X   = propagate_action(X,a(1,1),a(1,2),d,@(X,x,y)action_shift(X,x,y,d)) ;
        
%         fprintf('\nshift  %3d %3d %d : %d cones',x,y,d,sum(done))

end

X = transitive_closures( X ) ;
check_X(LL,X)

end