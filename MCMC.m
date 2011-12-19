function cone_map = MCMC( cone_map , ID )

if nargin>1 ,   cone_map.ID = ID ;     end

cone_map
cone_map.code.string    = file2str('MCMC.m') ;

default( cone_map , 'N_iterations'  , 100000)
default( cone_map , 'plot_every'    , 0     )
default( cone_map , 'plot_skip'     , 100   )
default( cone_map , 'display_every' , 50    )
default( cone_map , 'save_every'    , 200   )
default( cone_map , 'ID'            , 0     )
default( cone_map , 'max_time'      , 2000  )

% Initialize figure
if plot_every
scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)*0.7*0.5 1500*0.5 1200*1]) ;
end

% MAIN MCMC LOOP
fprintf('\n\nMCMC progress:')
t = cputime ;
tic

n_runs = 1 ;

X           = cone_map.initX ;
runbest     = X ;
runbest.i   = 1 ;
jj          = 1 ;
cone_map.bestX = {} ;
n_best      = 1 ;
while 1
    
    X = flip_MCMC( X, move(X, cone_map, [1 1]), @update_X, @(trial)trial.ll ) ;

    if X.ll>runbest.ll
        runbest = X ;
        runbest.i = jj ;
    end
        
    % reinitialize if stuck
    if jj - runbest.i > 400
        runbest = rmfield(runbest,{'invWW','masks','contact'}) ;
        cone_map.bestX{n_best} = runbest ;
        n_best  = n_best + 1 ;
        X       = cone_map.initX ;
        runbest = cone_map.initX ;
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

    if ~mod(jj,save_every) || jj>N_iterations || cputime-t>max_time
        cone_map.X          = X ;
        save(sprintf('result_%d',ID), 'cone_map' )
        if jj>N_iterations || cputime-t>max_time, break ; end
    end 
    jj = jj + 1 ;
end
fprintf('\n\ndone in %.1f sec\n\n',cputime - t) ;

end