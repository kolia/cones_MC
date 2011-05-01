function [bestX,cone_map] = MCMC( cone_map , ID )

if nargin>1 ,   cone_map.ID = ID ;     end

cone_map

default( cone_map , 'N_iterations'  , 100000)
default( cone_map , 'max_time'      , 2000  )
default( cone_map , 'plot_every'    , 0     )
default( cone_map , 'plot_skip'     , 100   )
default( cone_map , 'display_every' , 100   )
default( cone_map , 'save_every'    , 200   )
default( cone_map , 'ID'            , 0     )

% Initialize figure
if plot_every
scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)*0.7*0.5 1500*0.5 1200*0.5]) ;
end

% MAIN MCMC LOOP
fprintf('\n\nMCMC at temperature 1 :\n')
t = cputime ;
tic

n_runs = 1 ;

X           = cone_map.initX ;
runbest     = X ;
runbest.i   = 1 ;
jj = 1 ;
while 1
    
    X = flip_MCMC( move(X, cone_map, [1 1]), @update_X, @(trial)trial.ll ) ;

    if X.ll>runbest.ll
        runbest = X ;
        runbest.i = jj ;
    end
        
    % reinitialize if stuck
    if jj - runbest.i > 400
        bestX{n_runs} = runbest ;
        X       = initX ;
        runbest = initX ;
        runbest.i = jj ;
        n_runs  = n_runs + 1 ;
    end
    
    % DISPLAY stdout
    if ~mod(jj,display_every)
        fprintf('Iter%4d of %d  %4d cones    %6.0f L   %6.0f best   %8.2f sec\n',...
                            jj,N_iterations,numel(find(X.state>0)),X.ll,...
                            runbest.ll,toc)
        tic
    end
    
    % DISPLAY plot
    if ~mod(jj,plot_every)
        figure(h)
        plot_cones( X.state , cone_map ) ;
        title( sprintf('After %d MCMC iterations',jj),'FontSize' , 24 )
        % set(get(gca,'Title'),'Visible','on')
        drawnow
    end
    
    if ~mod(jj,save_every)
        bestX{n_runs} = runbest ;        
        save(sprintf('bestX_%d',ID), 'bestX')
    end
    
    jj = jj + 1 ;
    
    if jj>N_iterations || cputime-t>max_time ,  break ;  end
end
fprintf('\ndone in %.1f sec\n\n', cputime - t) ;

cone_map.X              = X ;
cone_map.bestX          = bestX ;
cone_map.code.MCMC      = file2str('MCMC.m') ;
save(sprintf('bestX_%d',ID), 'bestX')

end