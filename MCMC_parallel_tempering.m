function [cone_map,results,swap_stats] = MCMC_parallel_tempering( cone_map , ID )

if nargin>1
    cone_map.ID     = ID ;
end

betas           = cone_map.betas ;
N_iterations    = cone_map.N_iterations ;
N_instances     = length(betas) ;

default( cone_map , 'deltas' , ones(1,N_instances)) ;
default( cone_map , 'start_swap' , 0 )
default( cone_map , 'X'     , cell( N_instances , 1 ) )
default( cone_map , 'moves' , ...
                       [num2cell(1:N_instances) ...
                        num2cell(1:N_instances) ...
                        num2cell([1:N_instances-1 ; 2:N_instances],1)] )
default( cone_map , 'q' , 0.99 )

default( cone_map , 'plot_every'    , 0    )
default( cone_map , 'plot_skip'     , 100  )
default( cone_map , 'display_every' , 10   )
default( cone_map , 'save_every'    , 200  )
default( cone_map , 'track_instance', [1]  )
default( cone_map , 'ID'            , 0    )
default( cone_map , 'N_best'        , 3    )
default( cone_map , 'max_time'      , 1000  )


M0          = cone_map.M0 ;
M1          = cone_map.M1 ;
cone_map.SS = cone_map.cone_params.supersample ;
SS          = cone_map.SS ;
N_colors    = cone_map.N_colors ;

% sparse int matrix, with number of out-of-border adjacencies
cone_map.outofbounds = sparse([],[],[],M0*SS,M1*SS,2*(M0+M1)*SS) ;
cone_map.outofbounds(:,[1 M1*SS]) = 1 ;
cone_map.outofbounds([1 M0*SS],:) = cone_map.outofbounds([1 M0*SS],:) + 1 ;

fprintf('\n\nSTARTING %d MCMC instances with',N_instances)
fprintf('\n\n+ different inverse temperatures beta:\n')
fprintf('%4.2f  ',betas)
fprintf('\n\n+ different powers delta:\n')
fprintf('%4.2f  ',deltas)
fprintf('\n\n+ start swapping after:    %d iterations',start_swap)


% initializing variables
swap_stats      = cell( N_instances , 1 ) ;
results         = cell( N_instances , 1 ) ;
n_cones         = 0 ;

for i=1:N_instances
    
    % initialize X{i}
    if length(X)<i || isempty(X{i})
        X{i}            = initialize_X(M0,M1,N_colors,SS,...
                                cone_map.cone_params.replusion_radii,...
                                betas(i),deltas(i)) ;
    end
    
    X{i}.i          = i ;
    X{i}.SS         = cone_map.cone_params.supersample ;
    
    if isfinite( start_swap )
        swap_stats{i}.N50      = 0 ;
        swap_stats{i}.accepted = 0 ;
        swap_stats{i}.dSwap    = logical( sparse([],[],[],N_iterations,1) ) ;
    end
    
    if ismember(i,track_instance)
        results{i}.i         = i ;
        results{i}.iteration = 0 ;
        results{i}.swap      = logical( sparse([],[],[],N_iterations,1) ) ;
        results{i}.dX        = sparse([],[],[],N_iterations,3*X{1}.maxcones) ;
    end
    
end

bestX = cell(N_best,1) ;
for i=1:N_best
    bestX{i} = X{i} ;
end

if plot_every
scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)*0.7*0.5 1500*0.5 1200*0.5]) ;
end


% MAIN MCMC LOOP
fprintf('\n\nMCMC progress:\n')
t = cputime ;
tic

jj = 0 ;
while 1

    isswap = jj>start_swap ;

    for j=1:length(moves)
        this_move   = moves{j} ;
        i           = this_move(1) ;

        if isnumeric(this_move)

            % swap move if this_move has 2 indices
            if length(this_move) == 2 && isswap
                jjj = this_move(2) ;
                
%                 if isfield(swap_stats{i},'trials')          ...
%                 && X{i}.version == swap_stats{i}.version(1) ...
%                 && X{jjj}.version == swap_stats{i}.version(2) 
%                     swapX = swap_stats{i}.trials ;
%                     swap_stats{i} = rmfield(swap_stats{i},'trials') ;
% %                     fprintf('-')
%                 else
                    swapX = swap_closure( X{i} , X{jjj} , cone_map) ;
%                     fprintf('+')
%                 end
                
                swap_stats{i}.results{1} = results{i  } ;
                swap_stats{i}.results{2} = results{jjj} ;

                [ swap_stats{i} , swapX ] = flip_MCMC( ...
                    swap_stats{i}, swapX{1}, swapX(2:end), @update_swap ) ;
                
                results{i  } = swap_stats{i}.results{1} ;
                results{jjj} = swap_stats{i}.results{2} ;
                
                X{i} = swapX.X ; X{jjj} = swapX.with ;

            % regular MCMC move if this_move has one index
            elseif length(this_move) == 1
                [ results{i}, X{i} ] = flip_MCMC( ...
                    results{i}, X{i}, ...
                    move( X{i} , 2 , q , cone_map ), @update_X ) ;
            
            end
            n_cones = numel(find(X{1}.state>0)) ;
        end
    end

    if N_best 
        for inst=1:N_instances
            for b=1:N_best
                if X{inst}.ll>bestX{b}.ll && ...
                    ( b==1 || X{inst}.ll<bestX{b-1}.ll)
                    for bb=b+1:N_best
                        bestX{bb} = bestX{bb-1} ;
                    end
                    bestX{b} = X{inst} ;
                end
            end
        end
    end
    
    % DISPLAY stdout
    if ~mod(jj,display_every)
        if isswap
            swapstring = '   swaps' ;
        else
            swapstring = 'no swaps' ;
        end
        fprintf('Iteration:%4d of %d  %s  %4d cones %8.2f sec\n',...
                            jj,N_iterations,swapstring,n_cones,toc)
        if isswap
            fprintf('swap %% ')
            for ij=1:N_instances-1
                if ~isempty(swap_stats{ij})
                    fprintf('%3.0f ',...
                        100*swap_stats{ij}.accepted/swap_stats{ij}.N50)
%                         100*full( sum(results{ij}.dX(results{ij}.swap,1)>0)/...
%                               sum(results{ij}.swap) ))
                end
            end
            fprintf('   ')
        end
        tic
    end
    
    % DISPLAY plot
    if ~mod(jj,plot_every)
        figure(h)
        for i=1:plot_skip:max(1,(N_instances-plot_skip))
            NN = ceil(sqrt(N_instances/plot_skip)) ;
            subplot(NN,ceil(N_instances/(plot_skip*NN)),...
                                1+floor((i-1)/plot_skip))
            % subplot(1,N_instances,i)
            plot_cones( X{i}.state , cone_map ) ;
            titl = sprintf('\\beta %.2g',betas(i)) ;
            if i < N_instances && isswap
                s = swap_stats{i} ;
                titl = sprintf('%s     acc %.2f',titl,s.accepted/s.N50) ;
            end
            if i == 1 %ceil(N_instances/8)
                title({sprintf('Iteration %d',jj) ; titl },'FontSize',16)
            else
                title( titl , 'FontSize',16)
            end
            % set(get(gca,'Title'),'Visible','on')
        end
        drawnow
    end
    
    if ~mod(jj,save_every)
        allX = X ;
        X    = X(track_instance) ;
        save(sprintf('stats_%d',ID), 'results', 'swap_stats', 'X', 'bestX')
        X    = allX ;
        clear allX
    end
    
    if jj>N_iterations || cputime-t>max_time ,  break ;  end
end
fprintf('\ndone in %.1f sec\n\n',cputime - t) ;

cone_map.X              = X ;
cone_map.bestX          = bestX ;
cone_map.code.run_MCMC  = file2str('run_MCMC.m') ;
X           = X(track_instance) ;
save(sprintf('stats_%d',ID), 'results' , 'swap_stats', 'X', 'bestX')

end


function default( s , name , value )

if isfield( s , name )
    assignin( 'caller' , name , s.(name) ) ;
else
    assignin( 'caller' , name , value ) ;
end

end