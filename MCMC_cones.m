function cone_map = MCMC_cones( GC_stas , cone_params , cone_map )
%% cone_map = map_cones( GC_stas , cone_params [ , cone_map] )
%  Run MCMC to find cone locations.
%  All important parameters are within lines 17-26

if nargin<3  ,  cone_map = struct ;   end

% load variables defined in setup_cone_LL (shared with greedy_cones.m)
[STA_W,cone_map] = setup_cone_LL(GC_stas , cone_params , cone_map) ;
vars = fields(cone_map) ;
for i=1:length(vars)
    eval(sprintf('%s = cone_map.(vars{i}) ;',vars{i})) ;
end

%% prior_LL      prior log-likelihood : ~hard repulsion   VERY SENSITIVE
%   prior_LL      = @(X)-1e5*sum(sum( triu(X.overlaps,1) > max_overlap*0.01)) ; %   weak
  prior_LL      = @(X)-1e6*sum(sum( triu(X.overlaps,1) > max_overlap*0.01 )) ;  %   strong


%% APPROXIMATION!  sparsify STA_W
minSTAW = min(STA_W(:)) ;
STA_W(abs(STA_W)<minSTAW * 0.01 ) = 0 ;
STA_W = sparse(STA_W') ;

%% plot and display info every ? MCMC iterations  (0 for never)
plot_every      = 0 ;
display_every   = 100 ;


%% PARAMETERS FOR MCMC
  TOTAL_trials  = 1   * M2 * SS ; % number of trials after burn-in = TOTAL_trials * n_trials ;
  burn_in       = 30  * M2 * SS ; % number of burn-in trials

% q             probability of trying to move an existing cone vs. placing
%               a new one.
  q             = 0.8 ;

% these params shouldn't need tweaking unless the problem setup changes
%
% betas         number of independent instances run simultaneously
  betas         = ones(1,2) ;
%
% moves         sequence of moves at each iteration, currently:
%               - a regular MC move for each instance


%% MCMC RUN
N_instances     = length(betas) ;
n_trials        = 10 ;      % only applies after burn-in
fprintf('\n\nSTARTING %d MCMC instances with different inverse temperatures beta:\n',N_instances)
fprintf('%.2f   ',betas)

% move proposal distribution: higher probabilities where STAs were larger
MPD = reshape(full(STA_W),[sizeROI N_colors N_GC]) ;
MPD = MPD ./ repmat( max( max( max(abs(MPD),[],1) , [] , 2) , [] , 3 ) , [sizeROI N_colors 1] ) ;
MPD = sum( max(abs(MPD),[],4) ,3) ;

% initializing variables
proposal_prob   = cumsum(MPD(:)) ;
proposal_prob   = proposal_prob ./ proposal_prob(end) ;
flip_LL         = cell( N_instances , 1 ) ;
accumulated     = cell( N_instances , 1 ) ;
accumulator     = cell( N_instances , 1 ) ;
jitter          = cell( N_instances , 1 ) ;
X               = cell( N_instances , 1 ) ;

for i=1:N_instances
    X{i}.state      = sparse([],[],[],1,NROI,ceil(NROI/(15*SS^2))) ;
    accumulated{i}  = zeros( 1 , NROI + 2 ) ;
    accumulator{i}  = @(y)[ones(size(y,1),1) sum(y,2) y] ;
    jitter{i}       = @(X,n_trial)jitter_color(X , n_trial , q , sizeROI(1) , sizeROI(2) , proposal_prob , N_colors) ;

    flip_LL{i}      = @(X,flips)flip_color_LL( ...
        X , flips , prior_LL , cell_consts , STA_W , coneConv , colorDot , sizeROI , betas(i)) ;
    
    X{i} = flip_LL{i}( X{i} , [] ) ;  % initialize X{i}
end

% MC move sequence
% temp0   = repmat(2:N_instances-1,N_instances-2,1) ;
% temp1   = temp0' ;
% choose  = logical( tril(ones(N_instances-2),-1) ) ;
% pairs   = [temp1(choose) temp0(choose)]' ;

moves   = num2cell( 1:N_instances ) ;

% moves   = [ num2cell(  repmat( 1:N_instances , 1 , N_instances-2 ) )                ...
%             num2cell( pairs , 1)                                                    ...
%             num2cell( [ N_instances-1:-1:2 ; N_instances*ones(1,N_instances-2)] ,1) ...
%             num2cell( [         ones(1,N_instances-2) ; N_instances-1:-1:2] , 1 )   ] ;

N_moves      = length(moves) ;
burn_in      = ceil( burn_in / (n_trials * N_moves) ) ;
N_iterations = burn_in + ceil( TOTAL_trials / (n_trials * N_moves) ) ;
n_cones      = zeros( N_iterations , 1 ) ;

if plot_every
scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)*0.7 1500 600]) ;
end

% MAIN MCMC LOOP
% fprintf('\n\n      MCMC progress:|0%%              100%%|\n ')
fprintf('\n\nMCMC progress:')
t = cputime ;
tic
for jj=1:N_iterations
%     if ~mod(jj,floor(N_iterations/20)) , fprintf('*') , end

    isswap = jj>burn_in*0.8 ;
%     isswap = jj>9 ;

    for j=1:N_moves
        this_move = moves{j} ;

        if isnumeric(this_move)
            i = this_move(1) ;
            
            % swap move if this_move has 2 indices
            if length(this_move) == 2 && isswap
                swapX = swapper( X{i} , X{this_move(2)} ) ;
                swap_flipper = @(swapX,flips) swap_LL(swapX,flips, flip_LL{i} , flip_LL{this_move(2)}) ;

                if ~isempty(swapX.flips)
                    for ii=1:2
                    [ accumulated{i} , swapX ] = ...
                        flip_MCMC( accumulated{i} , swapX , accumulator{i} , @(X)X.flips , swap_flipper ) ;
                    end
                    X{i} = swapX.X ; X{this_move(2)} = swapX.with ;
                end
                
            % regular MCMC move if this_move has one index
            elseif length(this_move) == 1
                i = this_move ;
                if jj>burn_in
                    jit = @(X)jitter{i}(X,n_trials) ;
                else
                    jit = @(X)jitter{i}(X,1) ;
                end
                [ accumulated{i} , X{i} ] = ...
                    flip_MCMC( accumulated{i} , X{i} , accumulator{i} , jit , flip_LL{i} , jj<burn_in ) ;
            end
            n_cones(jj) = sum(X{1}.state) ;

            if jj>1 && n_cones(jj-1)>n_cones(jj)
                fprintf('\nN_CONES DECREASED!!!\n')
            end
        end
    end

    if ~mod(jj,display_every)
        if jj>burn_in
            swapstring = 'average' ;
        else
            swapstring = 'burn-in' ;
        end
        fprintf('\nIteration:%4d of %d \t %s\t  %4d cones \t%8.2f sec',...
                            jj,N_iterations,swapstring,n_cones(jj),toc)
        tic
    end
    
    if ~mod(jj,plot_every)
        figure(h)
        for i=1:N_instances
%             subplot(2,ceil(N_instances/2),i)
            subplot(1,N_instances,i)
            colormap('pink')
            GGG = zeros(26*SS,46*SS,3) ;
            GGG(ROI) = X{i}.state ;
            imagesc(GGG) ;
            titl = sprintf('X_%d   \\beta %.2f',i,betas(i)) ;
            if i == ceil(N_instances/4)
                title({sprintf('Iteration %d',jj) ; titl },'FontSize',22)
            else
                title( titl , 'FontSize',22)
            end
            %             set(get(gca,'Title'),'Visible','on')
        end
        drawnow
    end
    
end
fprintf('\ndone in %.1f sec\n\n',cputime - t) ;
cone_map.stas           = GC_stas ;
cone_map.cone_RF        = cone_RF ;
cone_map.cone_params    = cone_params ;
cone_map.prior_LL       = func2str(prior_LL) ;
cone_map.max_overlap    = max_overlap ;
cone_map.coneConv       = coneConv ;
cone_map.colorDot       = colorDot ;

ROI = logical(ROI) ;

% try
% for i=1:numel(cone_map.X)
%     best = cone_map.X{i}.best ;
%     for j=1:size(cone_map.X{1}.best,1)
%         cone_map.X{i}.best{j}.ll    = best(i,1) ;
%         cone_map.X{i}.best{j}.state = zeros(N,1) ;
%         cone_map.X{i}.best{j}.state(ROI) = best(i,2:end) ;
%     end
% end
% end

for i=1:N_instances
    cone_map.X{i}       = X{i} ;
    cone_map.X{i}.state = zeros(N,1) ;
    cone_map.X{i}.state(ROI) = X{i}.state ;
end
cone_map.betas          = betas ;
cone_map.n_cones        = n_cones ;
cone_map.burn_in        = burn_in ;
cone_map.N_iterations   = N_iterations ;
cone_map.moves          = moves ;
for i=1:N_instances
    cone_map.accumulated{i}              = zeros(N+2,1) ;
    cone_map.accumulated{i}(1:2)         = accumulated{i}(1:2) ;
    cone_map.accumulated{i}(find(ROI)+2) = accumulated{i}(3:end) ;
end
cone_map.ROI_super      = ROI ;
end