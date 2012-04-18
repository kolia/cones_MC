function to_save = CAST( cone_map , ID )

if nargin>1 ,   cone_map.ID = ID ;     end

cone_map.has_evidence = logical(squeeze(sum(abs(cone_map.LL),3))>0) ;
cone_map = rmfield(cone_map,{'LL','NICE'})
cone_map.code.string    = file2str('CAST.m') ;

default( cone_map , 'N_iterations'  , 1000000)
default( cone_map , 'plot_every'    , 0      )
default( cone_map , 'plot_skip'     , 100    )
default( cone_map , 'display_every' , 50     )
default( cone_map , 'save_every'    , 5000   )
default( cone_map , 'profile_every' , 0      )
default( cone_map , 'ID'            , 0      )
default( cone_map , 'max_time'      , 200000 )
default( cone_map , 'deltas' , ones(1,length(cone_map.betas))) ;
default( cone_map , 'N_fast'        , 1      )

% Initialize figure
if plot_every
scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)*0.7*0.5 1500*0.5 1200*1]) ;
end

fprintf('\n\nSTARTING CAST with' )
fprintf('\n\n+ different inverse temperatures beta:\n')
fprintf('%4.2f  ',cone_map.betas )
fprintf('\n\n+ different powers delta:\n')
fprintf('%4.2f  ',cone_map.deltas)

% initializing variables
X = cell(2,1) ;
for i=1:1+cone_map.N_fast ,  X{i} = remake_X(cone_map,cone_map.initX) ; end

X{1}.swap = logical( sparse([],[],[],N_iterations,1) ) ;
X{1}.dX   = sparse([],[],[],N_iterations,3*X{1}.maxcones) ;

% Tell update_X.m to save X{1}.ll whenever it updates X after an MCMC move.
X{1}.LL_history  = zeros(cone_map.N_iterations,1) ;

% Initialize array to save complete history ST.i in.
for j=1:cone_map.N_fast
    X{1+j}.STi_history = zeros(cone_map.N_iterations,1) ;
end

N_temp = length(cone_map.betas) ;
ST.T = cell(N_temp,1) ;
ST.i = ones(cone_map.N_fast,1) ;
ST.n = 1 ;
ST.g = exp(-3.1035+0.2268*(1:N_temp)) ; % from converged g of previous runs
for i=1:N_temp
    ST.T{i} = [cone_map.betas(i) cone_map.deltas(i)] ;
end

bestX     = X{1} ;

% MAIN MCMC LOOP
fprintf('\n\nMCMC progress:')
t = cputime ;
tic

% LL = @(x)get_LL(x,cone_map,T) ;

jj = 1 ;
while 1

    % regular MCMC for all instances
    X{1} = flip_MCMC( X{1}, move( X{1}, cone_map , [1 1]), cone_map, {[1 1]} ) ;

    for j=1:cone_map.N_fast
        X{1+j} = flip_MCMC( X{1+j}, move( X{1+j}, cone_map , ST.T{ST.i(j)}), cone_map, {ST.T{ST.i(j)}} ) ;
    end

    for j=1:cone_map.N_fast
        % swap move if X{1+j} is at T=1
        if ST.i(j)==1 && X{1}.N_cones>10  && X{1+j}.N_cones>10
            old_ll      = X{1}.ll ;
            [X{1},X{1+j}] = swap_step( X{1}, [1 1], X{1+j}, [1 1], cone_map ) ;
            if old_ll ~= X{1}.ll , fprintf(' swap dll : %.2f',X{1}.ll-old_ll) ; end
        end
        
        % update temperature step
        ST = SimTempMCMC( X{1+j}, cone_map, @get_LL, ST , j ) ;

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
        tic
    end

    % DISPLAY plot
    if ~mod(jj,plot_every)
        if ishandle(h)
            clf
            iters = max(1,X{1}.iteration-2000)+(1:2000) ;
            subplot(5,1,1:3)
            plot_cones( X{1}.state , cone_map ) ;
            title( sprintf('After %d MCMC iterations',jj),'FontSize' , 21 )
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
            ylabel( sprintf('fast Temps'),'FontSize' , 18 )
            xlabel('Iterations','FontSize' , 18)
            drawnow ; hold off
        end
    end
    
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
    
    if ~mod(jj,save_every) || jj>N_iterations || cputime-t>max_time
        to_save = rmfield(cone_map,{'STA','initX','sparse_struct'}) ;
%         cone_map.X          = X{1} ;
        to_save.X = rmfield(X{1},{'contact'}) ;
        try to_save.X = rmfield(X{1},{'invWW'}) ; end
%         try to_save.X = rmfield(X{1},{'WW'})    ; end
%         cone_map.bestX      = bestX ;
        to_save.bestX = rmfield(bestX,{'contact'}) ;
        try to_save.bestX = rmfield(bestX,{'invWW','LL_history','cputime','N_cones_history','dX','excluded','sparse_STA_W_state','swap'}) ; end
%         try to_save.bestX = rmfield(bestX,{'WW'})    ; end
        to_save.ST         = ST ;
        save(sprintf('result_%d',ID), 'to_save' )
        if jj>N_iterations || cputime-t>max_time, break ; 
        else clear to_save
        end
    end 
    jj = jj + 1 ;
end
fprintf('\n\ndone in %.1f sec\n\n',cputime - t) ;

end
