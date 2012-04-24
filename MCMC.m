function to_save = MCMC( cone_map , ID )

% defaults
default( cone_map , 'N_iterations'  , 1000000)
default( cone_map , 'plot_every'    , 0      )
default( cone_map , 'display_every' , 50     )
default( cone_map , 'save_every'    , 0      )
default( cone_map , 'ID'            , 0      )
default( cone_map , 'max_time'      , 200000 )
default( cone_map , 'save_disk_space', false  )

% IDs for each chain on the cluster; not useful for single local execution
if nargin>1 ,   cone_map.ID = ID ;     end

% has_evidence is true only where cones contribute to the likelihood
cone_map.has_evidence = logical(squeeze(sum(abs(cone_map.LL),3))>0) ;

% archive this file into saved result, for future reference
cone_map.code.string    = file2str('MCMC.m') ;

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

%% initialization with  'hot' greedy configuration
cone_map.track_contacts = true ;
greed_hot = GREEDY(cone_map, 'hot') ;
cone_map.initX = greed_hot.X ;

% % reduce memory footprint: LL is only used by greedy
% cone_map = rmfield(cone_map,'LL') ;

% initialize MCMC loop
X           = cone_map.initX ;
runbest     = X ;
runbest.i   = 1 ;
jj          = 1 ;
cone_map.bestX = {} ;
n_best      = 1 ;

% tell update_X to record accepted moves; used to replay chain and plots
X.dX   = sparse([],[],[],N_iterations,3*X.maxcones) ;

%% main MCMC loop
while 1
    
    % trials{i} is a proposed new configuration; each contains a new ll
    trials = move(X, cone_map, [1 1]) ;
    
    % choose one of the candidates, update_X its data structures
    X = flip_MCMC( X, trials, cone_map, {[1 1]} ) ;

    % keep track of best configuration encountered
    if X.ll>runbest.ll
        runbest = X ;
        runbest.i = jj ;
    end
        
    % reinitialize to cone_map.initX if MCMC becomes stuck
    if jj - runbest.i > cone_map.M0*cone_map.M1/3
        % use less disk space and memory: remove large data structures
        runbest = rmfield(runbest,{'masks','contact'}) ;
        try runbest = rmfield(runbest,'invWW') ; end
        try runbest = rmfield(runbest,{'LL_history','cputime','N_cones_history','dX','excluded','sparse_STA_W_state'}) ; end
        
        % record current best confguration before reinitializing
        cone_map.bestX{n_best} = runbest ;
        n_best  = n_best + 1 ;
        
        % reinitialize
        fprintf('reinitializing...\n')
        runbest = cone_map.initX ;
        runbest.LL_history = X.LL_history ;
        runbest.N_cones_history = X.N_cones_history ;
        runbest.cputime = X.cputime ;
        runbest.iteration = X.iteration ;
        runbest.i = jj ;
        n_runs  = n_runs + 1 ;
        X = runbest ;
    end
    
    % DISPLAY to stdout
    if ~mod(jj,display_every)
        fprintf('Iter%4d of %d  %4d cones    %6.0f L   %6.0f best   %8.2f sec\n',...
                            jj,N_iterations,numel(find(X.state>0)),X.ll,...
                            runbest.ll,toc)
        tic
    end

    % PLOT
    if ~mod(jj,plot_every)
        figure(h)
        plot_cones( X.state , cone_map ) ;
        title( sprintf('After %d MCMC iterations',jj),'FontSize' , 24 )
        % set(get(gca,'Title'),'Visible','on')
        drawnow
    end
    
    % SAVE to disk
    if ~mod(jj,save_every) || jj>N_iterations || cputime-t>max_time
        if save_disk_space
            % use less disk space: remove large data structures
            to_save = rmfield(cone_map,{'STA','initX','sparse_struct'}) ;
            to_save.X = rmfield(X,{'contact'}) ;
            try to_save.X = rmfield(X,{'invWW'}) ; end
        else
            to_save   = cone_map ;
            to_save.X = X ;
        end
        if ~mod(jj,save_every)
            save(sprintf('result_%d',ID), 'to_save' )
        end
        if jj>N_iterations || cputime-t>max_time, break ; 
        else clear to_save
        end
    end 
    jj = jj + 1 ;
end
fprintf('\n\ndone in %.1f sec\n\n',cputime - t) ;

end
