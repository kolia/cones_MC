function cone_map = MCMC_parallel_tempering( cone_map , ID )

if nargin>1 ,   cone_map.ID = ID ;     end

cone_map

default( cone_map , 'N_iterations'  , 100000)
default( cone_map , 'max_time'      , 2000  )
default( cone_map , 'plot_every'    , 0     )
default( cone_map , 'plot_skip'     , 100   )
default( cone_map , 'display_every' , 25    )
default( cone_map , 'save_every'    , 200   )
default( cone_map , 'ID'            , 0     )

default( cone_map , 'track_instance', [1]   )

% Initialize figure
if plot_every
scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)*0.7*0.5 1500*0.5 1200*0.5]) ;
end

N_iterations    = cone_map.N_iterations ;
N_instances     = length(cone_map.betas) ;

default( cone_map , 'deltas' , ones(1,N_instances)) ;
default( cone_map , 'start_swap' , 0 )
default( cone_map , 'X'     , cell( N_instances , 1 ) )
default( cone_map , 'moves' , ...
                       [num2cell(1:N_instances) ...
                        num2cell(1:N_instances) ...
                        num2cell([1:N_instances-1 ; 2:N_instances],1)] )

fprintf('\n\nSTARTING %d MCMC instances with',N_instances)
fprintf('\n\n+ different inverse temperatures beta:\n')
fprintf('%4.2f  ',cone_map.betas)
fprintf('\n\n+ different powers delta:\n')
fprintf('%4.2f  ',cone_map.deltas)
fprintf('\n\n+ start swapping after:    %d iterations',cone_map.start_swap)


for i=1:N_instances
    
    % initialize X{i}
    if length(X)<i || isempty(X{i})
        X{i} = cone_map.initX ;
        X{i}.T = [cone_map.betas(i) cone_map.deltas{i}] ;
    end    
    X{i}.i          = i ;
    X{i}.SS         = cone_map.cone_params.supersample ;
    
    if ismember(i,track_instance)
        X{i}.iteration = 0 ;
        X{i}.swap      = logical( sparse([],[],[],N_iterations,1) ) ;
        X{i}.dX        = sparse([],[],[],N_iterations,3*X{1}.maxcones) ;
    end
    
end

bestX     = X{1} ;
cone_map.initial_X = initial_X ;

% MAIN MCMC LOOP
fprintf('\n\nMCMC progress:\n')
t = cputime ;
tic

jj = 1 ;
while 1
    
    isswap = jj>start_swap ;

    for j=1:length(cone_map.moves)
        this_move   = cone_map.moves{j} ;
        i           = this_move(1) ;

        if isnumeric(this_move)

            ttt = cputime ;
            
            % swap move if this_move has 2 indices
            if length(this_move) == 2 && isswap
                jjj = this_move(2) ;

                [X{i},X{jjj}] = swap_step( X{i  } , [X{i  }.T],...
                                           X{jjj} , [X{jjj}.T],...
                                           cone_map) ;

            % regular MCMC move if this_move has one index
            elseif length(this_move) == 1
                X{i} = flip_MCMC( X{i}, move(X{i}, cone_map, X{i}.T), @update_X ) ;
            end
        end
        
    end

    if X{1}.ll>bestX.ll ,  bestX = X{1} ; end
    
    % DISPLAY stdout
    if ~mod(jj,display_every)
        if isswap
            swapstring = '   swaps' ;
        else
            swapstring = 'no swaps' ;
        end
        fprintf('Iteration:%4d of %d  %s  %4d cones %8.2f sec\n',...
                 jj,N_iterations,swapstring,numel(find(X{1}.state>0)),toc)
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
            titl = sprintf('\\beta %.2g',cone_map.betas(i)) ;
            if i == 1 %ceil(N_instances/8)
                title({sprintf('Iteration %d',jj) ; titl },'FontSize',16)
            else
                title( titl , 'FontSize',16)
            end
            % set(get(gca,'Title'),'Visible','on')
        end
        drawnow
    end
    
    if ~mod(jj,save_every) , save_XbestX(X(track_instance),bestX) ; end
    
    jj = jj + 1 ;
    
    if jj>N_iterations || cputime-t>max_time ,  break ;  end
    
end
fprintf('\ndone in %.1f sec\n\n',cputime - t) ;

cone_map.X              = X ;
cone_map.bestX          = bestX ;
cone_map.code.string    = file2str('MCMC_parallel_tempering.m') ;
save_XbestX(X(track_instance),bestX) ;
end


function save_XbestX( X , bestX )
save(sprintf('stats_%d',ID), 'X', 'bestX')
end