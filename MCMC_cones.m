function cone_map = MCMC_cones( GC_stas , cone_params , cone_map )
%% cone_map = map_cones( GC_stas , cone_params )
%  Run MCMC to find cone locations.
%  Assumes existence of ../resource folder, where STA_W

addpath(genpath(pwd))

% size of region of interest
  [M0,M1,N_colors] = size(GC_stas(1).spatial) ;
  M2 = M0*M1 ;

% supersampling factor
  SS          = cone_params.supersample ;
  
%% PARAMETERS FOR MCMC
  TOTAL_trials  = 20 * M2 * SS ;
  burn_in       = 10 * M2 * SS ;

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


%% PROBLEM SETUP: ganglion cell STAs

% set up Region of Interest
if isfield(cone_map,'ROI')
    ROI  = zeros(M0*SS,M1*SS,N_colors) ;
    for c=1:N_colors
        ROI(:,:,c) = kron(cone_map.ROI(:,:,c),ones(SS,SS)) ;
    end
else
    ROI  = ones(M0*SS,M1*SS,N_colors) ;
end

% stereotyped cone receptive field
s = cone_params.sigma ;
r = cone_params.support_radius ;
cone_RF = exp(-0.5 * ((-SS*r:SS*r)/(SS*s)).^2)' * exp(-0.5 * ((-SS*r:SS*r)/(SS*s)).^2) ;
cone_RF = cone_RF / sum(cone_RF(:)) ; clear r s

N = M2*(SS^2)*N_colors ;

% Unpacking GC_stas into: STA, norms of STAs and N_spikes
N_GC = length(GC_stas) ;
STA_norm = zeros(N_GC,1) ;
N_spikes = zeros(N_GC,1) ;
STA      = zeros(N,N_GC) ;
for i=1:N_GC
    N_spikes(i) = length(GC_stas(i).spikes) ;
    temp = cell(N_colors,1) ;
    for c=1:N_colors
        temp{c} = kron(GC_stas(i).spatial(:,:,c),ones(SS,SS)) ;
        temp{c} = temp{c}(:) ;
    end
    STA(:,i)    = cell2mat(temp) ;
    STA_norm(i) = norm(STA(:,i)) ;
end

cell_consts = N_spikes ./ exp(STA_norm/2) * cone_params.stimulus_variance ;

%% SETUP for Log-LIKELIHOOD calculations

% file = sprintf('../resources/STA_W_super%d.mat',SS) ;

% try load(file) ; fprintf('\nloaded %s\n',file) 
% catch
    % W = matrix representing all M2 possible cone receptive fields
    STA_W = zeros(N,N_GC) ;
    fprintf('Calculating STA_W  --  counting up to %d:\n',M0*SS)
    for i=1:M0*SS
        fprintf('%d ',i)
        for j=1:M1*SS
            for c=1:N_colors
                BW      = zeros(M0*SS,M1*SS) ;
                BW(i,j) = 1 ;
                BW      = imfilter(BW,cone_RF) ;
                
                filter  = kron(cone_params.colors(c,:),BW(:)') ;
                
                STA_W((c-1)*M0*M1*SS^2 + (j-1)*M0*SS + i , :) = filter * STA ;
            end
        end
    end
    STA_W = STA_W .* repmat( (N_spikes ./ cell_consts)' ,N,1) ;
    
%     save(file,'STA_W')
% end

coneConv    = conv2(cone_RF,cone_RF) ;
colorDot    = cone_params.colors * cone_params.colors' ;
max_overlap = coneConv(ceil(size(coneConv,1)/2),ceil(size(coneConv,2)/2)) ;



%% MCMC RUN
N_instances     = length(betas) ;
n_trials        = 10 ;
fprintf('\n\nSTARTING %d MCMC instances with different inverse temperatures beta:\n',N_instances)
fprintf('%.2f   ',betas)

% initializing variables
flat_probs      = ROI / sum(ROI(:)) ;
flip_LL         = cell( N_instances , 1 ) ;
accumulated     = cell( N_instances , 1 ) ;
accumulator     = cell( N_instances , 1 ) ;
jitter          = cell( N_instances , 1 ) ;
X               = cell( N_instances , 1 ) ;
prior_LL        = @(X)-1e12*sum(sum( triu(X.overlaps,1) > max_overlap*0.1 )) ;

for i=1:N_instances
    X{i}.state      = zeros( 1 , N     ) ;
    accumulated{i}  = zeros( 1 , N + 2 ) ;
    accumulator{i}  = @(y)[ones(size(y,1),1) sum(y,2) y] ;
    jitter{i}       = @(X)jitter_color( X , n_trials , M0*SS , M1*SS , flat_probs , 3 ) ;

    flip_LL{i}      = @(X,flips)flip_color_LL( ...
        X , flips , prior_LL , cell_consts , STA_W' , coneConv , colorDot , [M0 M1]*SS , betas(i)) ;
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


% X{1}.state(randi(length(X{1}.state),length(X{1}.state)*2/3,1)) = 1 ;
% 
% jitter{1}       = @(X)jitter_color( X , 1 , M0*SS , M1*SS , flat_probs , 3 ) ;
% for i=1:10
% [ accumulated{1} , X{1} ] = ...
%     flip_MCMC( accumulated{1} , X{1} , accumulator{1} , jitter{1} , flip_LL{1} ) ;
% end

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

end
fprintf('    done in %.1f sec\n\n',toc) ;
cone_map.stas           = GC_stas ;
cone_map.cone_RF        = cone_RF ;
cone_map.cone_params    = cone_params ;
cone_map.prior_LL       = func2str(prior_LL) ;
cone_map.max_overlap    = max_overlap ;
cone_map.coneConv       = coneConv ;
cone_map.colorDot       = colorDot ;

try
for i=1:numel(cone_map.X)
    best = cone_map.X{i}.best ;
    for j=1:size(cone_map.X{1}.best,1)
        cone_map.X{i}.best{j}.ll    = best(i,1) ;
        cone_map.X{i}.best{j}.state = best(i,2:end) ;
    end
end
end

cone_map.X              = X ;
cone_map.betas          = betas ;
cone_map.n_cones        = n_cones ;
cone_map.burn_in        = burn_in ;
cone_map.N_iterations   = N_iterations ;
cone_map.moves          = moves ;
cone_map.accumulated    = accumulated ;
cone_map.ROI_super   = ROI ;
end