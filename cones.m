function [cone_map,results] = cones( cone_map , ID )
%% cone_map = cones( cone_map , ID )
%  Run MCMC to find cone locations.

% addpath('../libraries/swaps')

cone_map    = rmfield(cone_map,'ROI') ;

M0          = cone_map.M0 ;
M1          = cone_map.M1 ;
cone_map.SS = cone_map.cone_params.supersample ;
SS          = cone_map.SS ;
N_colors    = cone_map.N_colors ;

% sparse int matrix, with number of out-of-border adjacencies
cone_map.outofbounds = sparse([],[],[],M0*SS,M1*SS,2*(M0+M1)*SS) ;
cone_map.outofbounds(:,[1 M1*SS]) = 1 ;
cone_map.outofbounds([1 M0*SS],:) = cone_map.outofbounds([1 M0*SS],:) + 1 ;

if nargin<2 , ID = 0 ; end
cone_map.ID = ID ;

%% plot and display info every ? MCMC iterations  (0 for never)
plot_every      = 0 ;
plot_skip       = 100 ;
display_every   = 10 ;
save_every      = 500 ;
track_every     = 1000 ;

%% PARAMETERS FOR MCMC
  N_iterations  = 6   * M0 * M1 ; % number of trials
  start_swap    = 1/4 * M0 * M1 ; % iteration where swapping starts

% betas         temperatures of independent instances run simultaneously
  betas         = make_deltas(0.1,1,2,5000) ;

% deltas        temperatures of independent instances run simultaneously
  deltas        = make_deltas(0.3,1,3,length(betas)) ;
  
% q             probability of trying to move an existing cone vs. placing
%               a new one.
  q             = 0.99 ;

% moves         sequence of moves at each iteration, currently:
%               - a regular MC move for each instance

cone_map.q = q ;

%% MCMC RUN
N_instances     = length(betas) ;

moves = [num2cell(1:N_instances) ...
         num2cell(1:N_instances) ...
         num2cell([1:N_instances-1 ; 2:N_instances],1)] ;
% moves = num2cell(1:N_instances) ;  % no swaps for now

N_moves      = length(moves) ;


fprintf('\n\nSTARTING %d MCMC instances with',N_instances)
fprintf('\n\n+ different inverse temperatures beta:\n')
fprintf('%4.2f  ',betas)
fprintf('\n\n+ different powers delta:\n')
fprintf('%4.2f  ',deltas)
fprintf('\n\n+ start swapping after:    %d iterations',start_swap)


% initializing variables
jitter          = cell( N_instances , 1 ) ;
X               = cell( N_instances , 1 ) ;
swap_stats      = cell( N_instances , 1 ) ;
results         = cell( N_instances , 1 ) ;
n_cones         = 0 ;

for i=1:N_instances
    
    jitter{i}       = @(X)move(X  , 2 , q , cone_map) ;
    
    % initialize X{i}
    X{i}            = initialize_X(M0,M1,N_colors,SS,...
                        cone_map.cone_params.replusion_radii,...
                        betas(i),deltas(i)) ;
    
    X{i}.i          = i ;
    
	X{i}.SS         = cone_map.cone_params.supersample ;
    
%     swap_stats{i}.N500     = 0 ;
    swap_stats{i}.N50      = 0 ;
    swap_stats{i}.accepted = 0 ;
    swap_stats{i}.dSwap    = logical( sparse([],[],[],N_iterations,1) ) ;
    
    if ~mod(i-1,track_every) || i==N_instances
%         results{i}.summed  = zeros( 1 , M0*SS*M1*SS*N_colors ) ;
%         results{i}.change  = false( 1 , M0*SS*M1*SS*N_colors ) ;
%         results{i}.times   = 0 ;
%         results{i}.N_iter  = 0 ;
%         results{i}.N500    = 0 ;
%         results{i}.LL_mean = 0 ;
%         results{i}.LL_mean_square = 0 ;
        results{i}.i         = i ;
        results{i}.iteration = 0 ;
        results{i}.swap      = logical( sparse([],[],[],N_iterations,1) ) ;
        results{i}.dX        = sparse([],[],[],N_iterations,3*X{1}.maxcones) ;
    end
    
end


if plot_every
scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)*0.7*0.5 1500*0.5 1200*0.5]) ;
end

% hswaps = figure ;

GG0 = cone_map.NICE .^ 0.7 ; % !!! enhance low values visually


% MAIN MCMC LOOP
% fprintf('\n\n      MCMC progress:|0%%              100%%|\n ')
fprintf('\n\nMCMC progress:\n')
t = cputime ;
tic
for jj=1:N_iterations
%     if ~mod(jj,floor(N_iterations/20)) , fprintf('*') , end

    isswap = jj>start_swap ;

    for j=1:N_moves
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
                    results{i}, X{i}, jitter{i}(X{i}), @update_X) ;
            
            end
            n_cones = numel(find(X{1}.state>0)) ;
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
            colormap('pink')
            GGG = GG0 ;
            for c=1:3
                [ix,iy] = find(X{i}.state == c) ;
                for cc=1:3
                    GGG(ix + M0*SS*(iy-1) + M0*SS*M1*SS*(cc-1)) = 0 ;
                end
                GGG(ix + M0*SS*(iy-1) + M0*SS*M1*SS*(c-1)) = 1 ;
            end
            imagesc(GGG) ;
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

%         for i=1:N_instances
%             if isfield(results{i},'summed')
%                 cone_map.stats{i}   = rmfield(results{i},'change') ;
%                 cone_map.stats{i}.summed = ...
%                     reshape( results{i}.summed, [M0*SS M1*SS N_colors]) ;
%                 cone_map.stats{i}.N_iter = results{i}.N_iter ;
%             end
%         end
%         
%         MAX = max(cone_map.stats{1}.summed(:)) ;
%         if MAX>0
%             acc = cone_map.stats{1}.summed ./ MAX ;
%             imagesc(acc)
%             titl = sprintf('X_%d   \\beta %.2g',i,betas(1)) ;
%             if i == ceil(N_instances/4)
%                 title({sprintf('Iteration %d',jj) ; titl },'FontSize',16)
%             else
%                 title( titl , 'FontSize',16)
%             end
%             drawnow
%         end
    end
    
    if ~mod(jj,save_every)
        save(sprintf('stats_%d',ID), 'results' , 'swap_stats')
    end
end
fprintf('\ndone in %.1f sec\n\n',cputime - t) ;

save(sprintf('stats_%d',ID), 'results' , 'swap_stats', 'X')

cone_map.X              = X ;
cone_map.code           = file2str('cones.m') ;
cone_map.betas          = betas ;
cone_map.n_cones        = n_cones ;
cone_map.N_iterations   = N_iterations ;
cone_map.moves          = moves ;
% for i=1:N_instances
%     if isfield(results{i},'summed')
%         cone_map.stats{i}   = rmfield(results{i},'change') ;
%         cone_map.stats{i}.summed = ...
%             reshape( results{i}.summed , [M0*SS M1*SS N_colors] ) ;
%         cone_map.stats{i}.N_iter = results{i}.N_iter ;
%     end
% end
cone_map.swap_stats     = swap_stats ;

% if plot_every
%     figure
%     acc = cone_map.stats{1}.summed ./ max(cone_map.stats{1}.summed(:)) ;
%     imagesc(acc)
%     titl = sprintf('X_%d   \\beta %.2g',i,betas(1)) ;
%     if i == ceil(N_instances/4)
%         title({sprintf('Iteration %d',jj) ; titl },'FontSize',16)
%     else
%         title( titl , 'FontSize',16)
%     end
% end

end