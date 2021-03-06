function to_save = CAST( cone_map , ID )

% IDs for each chain on the cluster; not useful for single local execution
if nargin>1 ,   cone_map.ID = ID ;     end

% has_evidence is true only where cones contribute to the likelihood
cone_map.has_evidence = logical(squeeze(sum(abs(cone_map.LL),3))>0) ;

% archive this file into saved result, for future reference
cone_map.code.string    = file2str('CAST.m') ;

% defaults
default( cone_map , 'N_iterations'      , 1000000)
default( cone_map , 'max_time'          , 200000 )
default( cone_map , 'plot_every'        , 1000   )
default( cone_map , 'display_every'     , 50     )
default( cone_map , 'save_every'        , 0      )
default( cone_map , 'profile_every'     , 0      )
default( cone_map , 'ID'                , 0      )
default( cone_map , 'N_fast'            , 1      )
default( cone_map , 'swap_N_times'      , 50     )
default( cone_map , 'save_disk_space'   , false  )

% default hottest inverse temperature, max number of temps, and curvature
default( cone_map , 'min_beta'          , 0.2    )
default( cone_map , 'min_delta'         , 0.2    )
default( cone_map , 'max_temps'         , 130    )
default( cone_map , 'curvature'         , 0.1    ) % how do temps fall off?

% Initialize figure
if plot_every
scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)*0.7*0.5 1500*0.5 1200*1]) ;
end

%% initializations

% initialize with  'hot' greedy configuration
cone_map.track_contacts = true ;
greed_hot = GREEDY(cone_map, 'hot') ;
cone_map.initX = greed_hot.X ;

% % reduce memory footprint: LL is only used by greedy
% cone_map = rmfield(cone_map,'LL')

% initialize slow chain X{1} and fast chain X{2}
X = cell(1+cone_map.N_fast,1) ;
for i=1:1+cone_map.N_fast ,  X{i} = cone_map.initX ; end

initial_iterations = cone_map.initX.iteration ;

% initialize current best configuration
bestX = X{1} ;

% tell udate_swap to record iterations where a swap was accepted
X{1}.swap = logical( sparse([],[],[],N_iterations,1) ) ;

% tell update_X to record accepted moves; used to chain replay and plots
X{1}.dX   = sparse([],[],[],N_iterations,3*X{1}.maxcones) ;

% tell update_X to record X{1}.ll at each iteration
X{1}.LL_history  = zeros(cone_map.N_iterations,1) ;

% ST.T = the fixed progression of temperatures, ST.i = current temps
ST.T = cell(2,1) ;
ST.i = ones(cone_map.N_fast,1) ;
ST.n = 1 ;
ST.curvature = curvature ;
ST.T{1} = [1 1] ;
ST.T{2} = [min_beta min_delta] ;

% ST.g is used by Wang-Landau scheme in SimTempMCMC.m
ST.lg = ones(1,2) ;
ST.max_temps = max_temps ;

% X{2}.STi_history records the complete history of ST.i
for j=1:cone_map.N_fast
    X{1+j}.STi_history = zeros(cone_map.N_iterations,1) ;
end


%% MAIN CAST LOOP
fprintf('\nSTARTING CAST with' )
fprintf('\n+ hottest inverse temperature beta: %f\n',min_beta )
fprintf('+ smallest power delta            : %f\n',min_delta)
fprintf('\n\nCAST progress:')
t = cputime ;
tic

jj = 1 ;
while 1

    % regular MCMC move at temperature [1 1] for slow chain X{1}
    proposal = move( X{1}, cone_map , [1 1]) ;
    accept   = metropolis_hastings( X{1}.ll, proposal.ll, proposal.proposal_bias ) ;
    X{1} = update_X( {X{1}; proposal}, accept+1 ) ;

    % regular MCMC move at temperature ST.T{ST.i(j)} for fast chain X{2}
    for j=1:cone_map.N_fast
        proposal = move( X{1+j}, cone_map , ST.T{ST.i(j)}) ;
        accept   = metropolis_hastings( X{1}.ll, proposal.ll, proposal.proposal_bias ) ;
        X{1+j} = update_X( {X{1+j}; proposal}, accept+1 ) ;
    end
    
    for j=1:cone_map.N_fast
        % swap move if X{1+j} is at T=1
        if isfield(ST,'k') && 2*numel(ST.T)-1>max_temps && ...
                   ST.i(j)==1 && X{1}.N_cones>10  && X{1+j}.N_cones>10
            old_ll   = X{1}.ll ;
            old_fast = X{2}.ll ;
            [X{1},X{1+j}] = swap_step(X{1},X{1+j},cone_map, swap_N_times, ST.T{1}) ;
%             if old_ll ~= X{1}.ll
%                 fprintf(' dLL : %.2f',X{1}.ll-old_ll) ;
%             end
%             fprintf('  old_ll: %.2f  slow : %.2f  fast : %.2f  ',old_ll,X{1}.ll,X{2}.ll) ; 
            fprintf('  SWAPDLL: %.2f, SLOWDLL: %.2f',...
                X{1}.ll+X{2}.ll-old_ll-old_fast,...
                X{1}.ll-old_ll) ;
        end
        
        % update temperature step
        [X{1+j},ST] = SimTempMCMC( X{1+j}, cone_map, ST , j ) ;

        % save current ST.STi to X{2}.STi_history
        X{1+j}.STi_history(X{1+j}.iteration) = ST.i(j) ;
    end
     
    if X{1}.ll>bestX.ll ,  bestX = X{1} ; end

    % DISPLAY stdout
    if ~mod(jj,display_every)
        fprintf('\nIteration:%4d of %d  %4d cones %6.0f  ST.i: ', ...
                     jj,N_iterations,numel(find(X{1}.state>0)),X{1}.ll)
        fprintf('%2d ',ST.i)
        fprintf('%6.2f sec',toc)
        props = numel(ST.f)*ST.f/sum(ST.f) ;
%         fprintf('\nST.lg:') ; fprintf(' %g',ST.lg)
%         fprintf('\nST.f:') ;  fprintf(' %g',props)
        fprintf('   min ST.f %g >? 0.8',min(props))
        tic
    end

    % DISPLAY plot
    if ~mod(jj,plot_every)
        if ishandle(h)
            clf
            iters = max(initial_iterations,X{1}.iteration-2000)+(1:2000) ;
            iters = iters(X{1}.LL_history(iters)>0) ;
            subplot(5,1,1:3)
            plot_cones_matlab( X{1}.state , cone_map ) ;
            title( sprintf('After %d CAST iterations',jj),'FontSize' , 21 )
            subplot(5,1,4) ;
            plot(iters,X{1}.LL_history(iters))
            ylabel( sprintf('X(1) LL'),'FontSize' , 18 )
            subplot(5,1,5) ;
            hold on
            plotme = zeros(cone_map.N_fast,length(iters)) ;
            for j=1:cone_map.N_fast
                plotme(j,:) = X{1+j}.STi_history(iters) ;
            end
            plot( repmat(iters,cone_map.N_fast,1)' , plotme' )
            ylabel( sprintf('fast temp.'),'FontSize' , 18 )
            xlabel('Iterations','FontSize' , 18)
            drawnow ; hold off
        end
    end
    
    % PROFILING info saved to disk every profile_every iterations
    if profile_every
        if jj==1
            profile clear
            profile on
        elseif ~mod(jj,profile_every)
            p = profile('info');
            save(sprintf('profdat_%d',floor(jj/profile_every)),'p')
            profile clear
        end
    end
    
    % SAVE to disk
    if ~mod(jj,save_every) || jj>N_iterations || cputime-t>max_time
        if save_disk_space
            % use less disk space: remove large data structures
            to_save = rmfield(cone_map,{'STA','initX','sparse_struct'}) ;
            to_save.X = rmfield(X{1},{'contact'}) ;
            try to_save.X = rmfield(X{1},{'invWW'}) ; end
        else
            to_save = cone_map ;
            to_save.X = X ;
        end
        to_save.bestX = rmfield(bestX,{'contact'}) ;
        try to_save.bestX = rmfield(bestX,{'invWW','LL_history','cputime',...
            'N_cones_history','dX','excluded','sparse_STA_W_state','swap'}) ; end
        to_save.ST         = ST ;
        if ~mod(jj,save_every)
            save(sprintf('cast_result_%d',ID), 'to_save' )
        end
        if jj>N_iterations || cputime-t>max_time, break ; 
        else clear to_save
        end
    end 
    jj = jj + 1 ;
end
fprintf('\n\ndone in %.1f sec\n\n',cputime - t) ;

end
