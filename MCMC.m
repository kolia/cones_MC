function to_save = MCMC( cone_map , ID )

if nargin>1 ,   cone_map.ID = ID ;     end

cone_map.has_evidence = logical(squeeze(sum(abs(cone_map.LL),3))>0) ;
cone_map = rmfield(cone_map,{'LL','NICE'})
cone_map.code.string    = file2str('MCMC.m') ;

default( cone_map , 'N_iterations'  , 1000000)
default( cone_map , 'plot_every'    , 0      )
default( cone_map , 'plot_skip'     , 100    )
default( cone_map , 'display_every' , 50     )
default( cone_map , 'save_every'    , 5000   )
default( cone_map , 'ID'            , 0      )
default( cone_map , 'max_time'      , 200000 )

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

cone_map.initX = remake_X(cone_map,cone_map.initX) ;
X           = cone_map.initX ;
runbest     = X ;
runbest.i   = 1 ;
jj          = 1 ;
cone_map.bestX = {} ;
n_best      = 1 ;

while 1
    
    trials = move(X, cone_map, [1 1]) ;
    X = flip_MCMC( X, trials, cone_map, {[1 1]} ) ;
    
    if X.ll>runbest.ll
        runbest = X ;
        runbest.i = jj ;
    end
        
    % reinitialize if stuck
    if jj - runbest.i > cone_map.M0*cone_map.M1/3
        runbest = rmfield(runbest,{'masks','contact'}) ;
        try runbest = rmfield(runbest,{'invWW'}) ; end
        try runbest = rmfield(runbest,{'WW'})    ; end
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
%         cone_map.X          = X ;
        to_save = rmfield(cone_map,{'STA','initX','sparse_struct'}) ;
        to_save.X = rmfield(X,{'contact'}) ;
        try to_save.X = rmfield(X,{'invWW'}) ; end
%         try to_save.X = rmfield(X,{'WW'})    ; end
        save(sprintf('result_%d',ID), 'to_save' )
        if jj>N_iterations || cputime-t>max_time, break ; 
        else clear to_save
        end
    end 
    jj = jj + 1 ;
end
fprintf('\n\ndone in %.1f sec\n\n',cputime - t) ;

end
