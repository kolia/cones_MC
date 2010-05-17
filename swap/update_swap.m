function X = update_swap(X,trial,XLL,OLL)

% check_X(XLL,X.X)
% check_X(OLL,X.with)

% type = trial.move{1} ;
a    = trial.move{2} ;  Na = length(a) ;
b    = trial.move{3} ;  Nb = length(b) ;
LL_X = trial.move{4} ;
LL_O = trial.move{5} ;

% switch type    
%     case 'patch'
%         fprintf('\npatch  %3d %3d    N_cones %d %d',Na,Nb,X.X.N_cones,X.with.N_cones)
        
        ax = zeros(Na,1) ;
        ay = zeros(Na,1) ;
        ac = zeros(Na,1) ;
        for i=1:Na
            [x,y]   = find( X.X.id == a(i) , 1 ) ;
            ax(i)   = x ;
            ay(i)   = y ;
            ac(i)   = X.X.state(x,y) ;
            X.X     = update_cone_deletion(X.X,x,y) ;
        end
        
        bx = zeros(Nb,1) ;
        by = zeros(Nb,1) ;
        bc = zeros(Nb,1) ;
        for i=1:Nb
            [x,y]   = find( X.with.id == b(i) , 1 ) ;
            bx(i)   = x ;
            by(i)   = y ;
            bc(i)   = X.with.state(x,y) ;
            X.with  = update_cone_deletion(X.with,x,y) ;
        end
        
        for i=1:Na
            X.with  = update_cone_addition(X.with,ax(i),ay(i),ac(i),OLL) ;
        end
        
        for i=1:Nb
            X.X     = update_cone_addition(X.X   ,bx(i),by(i),bc(i),XLL) ;
        end
        
        X.X     = transitive_closures( X.X    ) ;
        X.with  = transitive_closures( X.with ) ;

        X.X.ll    = LL_X ;
        X.with.ll = LL_O ;
        
%     case 'whole'
%         fprintf('\nwhole  %3d %3d',Na,Nb)
% end

X.ll = trial.ll ;
% check_X(XLL,X.X)
% check_X(OLL,X.with)

end