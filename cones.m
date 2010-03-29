function cone_map = cones( cone_map )
%% cone_map = cones( LL )
%  Run MCMC to find cone locations.

% addpath('../library/kdtree')

LL = cone_map.LL ;
N_cones_factor = cone_map.N_cones_factor ;

[M0,M1,N_colors] = size(LL) ;

%% plot and display info every ? MCMC iterations  (0 for never)
plot_every      = 20 ;
display_every   = 20 ;


%% PARAMETERS FOR MCMC
  TOTAL_trials  = 5   * M0 * M1 ; % number of trials after burn-in = TOTAL_trials * n_trials ;
  burn_in       = 5   * M0 * M1 ; % number of burn-in trials

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

% initializing variables
flip_LL         = cell( N_instances , 1 ) ;
accumulated     = cell( N_instances , 1 ) ;
accumulator     = cell( N_instances , 1 ) ;
jitter          = cell( N_instances , 1 ) ;
X               = cell( N_instances , 1 ) ;

exclusion       = 10 ;

dprior          = @(X,x,y,c)dist_prior(X,x,y,c, exclusion, ...
                        @(d)-1e8*(abs(d-exclusion^2/2)<exclusion^2/2)  , ...
                        @(n) 0 ) ; % -10*(10-n) )  ;

for i=1:N_instances
    
    accumulated{i}  = zeros( 1 , M0*M1 + 1) ;
    accumulator{i}  = @(X)[1 X.state(:)'] ;
        
    flip_LL{i}      = @(X,x,y,c)flip_diag_LL( X , x , y , c , LL , N_cones_factor , dprior ) ;
    
    jitter{i}       = @(X,n_trial)jitter_color(X , n_trial , q , flip_LL{i} , N_colors) ;

    X{i}            = flip_LL{i}( struct , 0 , 0 , 0 ) ;  % initialize X{i}
    
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

GG0 = LL - min(LL(:)) ;
GG0 = GG0 / max(GG0(:)) ;

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
                    flip_MCMC( accumulated{i} , X{i} , accumulator{i} , jit , jj<burn_in ) ;
            end
            n_cones(jj) = numel(find(X{1}.state>0)) ;

%             if jj>1 && n_cones(jj-1)>n_cones(jj)
%                 fprintf('\nN_CONES DECREASED!!!\n')
%             end
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
            GGG = GG0 ;
            for c=1:3
                [ix,iy] = find(X{i}.state == c) ;
                for cc=1:3
                    GGG(ix + M0*(iy-1) + M0*M1*(cc-1)) = 0 ;
                end
                GGG(ix + M0*(iy-1) + M0*M1*(c-1)) = 1 ;
            end
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
cone_map.dprior       = func2str(dprior) ;

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
end
cone_map.betas          = betas ;
cone_map.n_cones        = n_cones ;
cone_map.burn_in        = burn_in ;
cone_map.N_iterations   = N_iterations ;
cone_map.moves          = moves ;
cone_map.accumulated    = accumulated ;

end