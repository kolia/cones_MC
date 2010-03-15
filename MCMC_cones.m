function cone_map = MCMC_cones( GC_stas , cone_params , cone_map )
%% cone_map = map_cones( GC_stas , cone_params )
%  Run MCMC to find cone locations.
%  Assumes existence of ../resource folder, where STA_W

if nargin<3  ,  cone_map = struct ;   end

% load variables defined in setup_cone_LL
[STA_W,cone_map] = setup_cone_LL(GC_stas , cone_params , cone_map) ;
vars = fields(cone_map) ;
for i=1:length(vars)
    eval(sprintf('%s = cone_map.(vars{i}) ;',vars{i})) ;
end


%% PARAMETERS FOR MCMC
  TOTAL_trials  = 200 * M2 * SS ;
  burn_in       = 100 * M2 * SS ;

% these params shouldn't need tweaking unless the problem setup changes
%        
% prior_LL      prior used both for greedy and MCMC
%
% betas         instance X_i has inverse temperature betas(i)
  betas         = [1 0.25 0.2 0.2 0.2 0.15] ;
%               X_1 is the main instance
%               X_(2..N-1) are a pool of warm instances
%               X_N is a hot instance that provides randomness to the pool
%
% moves         sequence of moves at each iteration, currently:
%               - 3 regular MC moves for each instance
%               - one swap between X_1 and each X_(2..N-1)
%               - one swap among each pair of   X_(2..N-1)
%               - one swap between X_N and each X_(2..N-1)


%% WHAT MCMC IS DOING: the short version 
% The MCMC used here is a combination of parallel tempering, aka replica 
% exchange MC, and cluster MC, with clusters adapted to our problem.
%
% 6 MC instances are run in parallel: one with the actual log-likelihood LL
% the other 5 instances i with log-likelihood LL * betas(i) < LL.
% Once in while, connected components of the symmetric difference between 
% pairs of instance configurations are swapped between instances.


%% MCMC RUN
N_instances     = length(betas) ;
n_trials        = 10 ;
fprintf('\n\nSTARTING %d MCMC instances with different inverse temperatures beta:\n',N_instances)
fprintf('%.2f   ',betas)

% initializing variables
flat_probs      = ones(NROI,1) / NROI ;
flip_LL         = cell( N_instances , 1 ) ;
accumulated     = cell( N_instances , 1 ) ;
accumulator     = cell( N_instances , 1 ) ;
jitter          = cell( N_instances , 1 ) ;
X               = cell( N_instances , 1 ) ;
prior_LL        = @(X)-1e12*sum(sum( triu(X.overlaps,1) > max_overlap*0.1 )) ;

for i=1:N_instances
    X{i}.state      = zeros( 1 , NROI     ) ;
    accumulated{i}  = zeros( 1 , NROI + 2 ) ;
    accumulator{i}  = @(y)[ones(size(y,1),1) sum(y,2) y] ;
    jitter{i}       = @(X)jitter_color( X , n_trials , sizeROI(1) , sizeROI(2) , flat_probs , 3 ) ;

    flip_LL{i}      = @(X,flips)flip_color_LL( ...
        X , flips , prior_LL , cell_consts , STA_W' , coneConv , colorDot , sizeROI , betas(i)) ;
end

% MC move sequence
temp0   = repmat(2:N_instances-1,N_instances-2,1) ;
temp1   = temp0' ;
choose  = logical( tril(ones(N_instances-2),-1) ) ;
pairs   = [temp1(choose) temp0(choose)]' ;

moves   = [ num2cell(  repmat( 1:N_instances , 1 , N_instances-2 ) )                ...
            num2cell( pairs , 1)                                                    ...
            num2cell( [ N_instances-1:-1:2 ; N_instances*ones(1,N_instances-2)] ,1) ...
            num2cell( [         ones(1,N_instances-2) ; N_instances-1:-1:2] , 1 )   ] ;

N_moves      = length(moves) ;
burn_in      = ceil( burn_in / (n_trials * N_moves) ) ;
N_iterations = burn_in + ceil( TOTAL_trials / (n_trials * N_moves) ) ;
n_cones      = zeros( N_iterations , 1 ) ;


% MAIN MCMC LOOP
fprintf('\n\n      MCMC progress\n|0%%              100%%|\n ')
tic
for jj=1:N_iterations
    if ~mod(jj,floor(N_iterations/20)) , fprintf('*') , end
    
    for j=1:N_moves
        this_move = moves{j} ;

        if isnumeric(this_move)
            i = this_move(1) ;
            
            % swap move if this_move has 2 indices
            if length(this_move) == 2 && jj>burn_in*0.8
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
                [ accumulated{i} , X{i} ] = ...
                    flip_MCMC( accumulated{i} , X{i} , accumulator{i} , jitter{i} , flip_LL{i} ) ;
            end
            n_cones(jj) = sum(X{1}.state) ;
        end
    end
    
%     GGG = zeros(26,46,3) ;
%     GGG(ROI) = X{1}.state ;
%     imagesc(GGG) ;
%     drawnow

end
fprintf('    done in %.1f sec\n\n',toc) ;
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