function cone_map = cones( cone_map )
%% cone_map = cones( LL )
%  Run MCMC to find cone locations.

% addpath('../library/kdtree')

LL = cone_map.LL ;
N_cones_factor = cone_map.N_cones_factor ;

[M0,M1,N_colors] = size(LL) ;

%% plot and display info every ? MCMC iterations  (0 for never)
plot_every      = 0 ;
display_every   = 50 ;


%% PARAMETERS FOR MCMC
  TOTAL_trials  = 100 * M0 * M1 ; % number of trials after burn-in = TOTAL_trials * n_trials ;
  burn_in       = 80  * M0 * M1 ; % number of burn-in trials

% q             probability of trying to move an existing cone vs. placing
%               a new one.
  q             = 0.8 ;

% these params shouldn't need tweaking unless the problem setup changes
%
% betas         number of independent instances run simultaneously
  betas         = [1 0.9 0.8 0.6 0.4 0.2] ;
%
% exclusions    exclusion distance between cones of instance i
  exclusions    = [9 8.98 8.96 8.85 8.7 8.58] ;
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
jitter          = cell( N_instances , 1 ) ;
X               = cell( N_instances , 1 ) ;
dprior          = cell( N_instances , 1 ) ;

N_cones_prior = @(n,i) -20*(i>1)*n*betas(i) ;

for i=1:N_instances
    

    dprior{i}       = @(X,x,y,c) betas(i) * dist_prior(X,x,y,c, exclusions(i), ...
                        @(d)-1e8*(abs(d-exclusions(i)^2/2)<exclusions(i)^2/2)  , ...
                        @(n)N_cones_prior(n,i) )  ;

    flip_LL{i}      = @(X,x,y,c)flip_diag_LL( X , x , y , c , LL * betas(i) , N_cones_factor * betas(i) , dprior{i} ) ;
    
    jitter{i}       = @(X,n_trial)jitter_color(X , n_trial , q , flip_LL{i} , N_colors) ;

    X{i}            = flip_LL{i}( struct , 0 , 0 , 0 ) ;  % initialize X{i}
    
end

accumulator  = @(X)[(X.state(:)'==1) (X.state(:)'==2) (X.state(:)'==3)] ;
accumulated  = zeros( 1 , M0*M1*N_colors ) ;

% MC move sequence
% temp0   = repmat(2:N_instances-1,N_instances-2,1) ;
% temp1   = temp0' ;
% choose  = logical( tril(ones(N_instances-2),-1) ) ;
% pairs   = [temp1(choose) temp0(choose)]' ;
% 
% moves   = [ num2cell(  repmat( 1:N_instances , 1 , N_instances-2 ) )                ...
%             num2cell( pairs , 1)                                                    ...
%             num2cell( [ N_instances-1:-1:2 ; N_instances*ones(1,N_instances-2)] ,1) ...
%             num2cell( [         ones(1,N_instances-2) ; N_instances-1:-1:2] , 1 )   ] ;

% moves   = num2cell( 1:N_instances ) ;

% moves   = [ num2cell(  1:N_instances )                          ...
%             num2cell( [1:N_instances-1 ; 2:N_instances] , 1 )   ] ;

% moves = {1 2 3 4 [1 2] [1 3] [1 4]} ;
% moves = {1 2 3 4 5 6 [1 2] [1 3] [2 3] [2 4] [3 4] [3 5] [4 5] [4 6] [5 6]} ;
moves = {1 2 3 4 5 6 [1 2] [2 3] [3 4] [4 5] [5 6]} ;

N_moves      = length(moves) ;
burn_in      = ceil( burn_in / (n_trials * N_moves) ) ;
N_iterations = burn_in + ceil( TOTAL_trials / (n_trials * N_moves) ) ;
n_cones      = zeros( N_iterations , 1 ) ;

stats = cell(N_moves,1) ;
for j=1:N_moves
    stats{j}.N        = 0 ;
    stats{j}.accepted = 0 ;
end

if plot_every
scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)*0.7 1500 1200]) ;
end

% hswaps = figure ;

GG0 = LL - min(LL(:)) ;
GG0 = ( GG0 / max(GG0(:))) .^ 0.8 ; % !!! take power to enhance low values visually

% MAIN MCMC LOOP
% fprintf('\n\n      MCMC progress:|0%%              100%%|\n ')
fprintf('\n\nMCMC progress:')
t = cputime ;
tic
for jj=1:N_iterations
%     if ~mod(jj,floor(N_iterations/20)) , fprintf('*') , end

    isswap = jj>burn_in*0.1 ;

    for j=1:N_moves
        this_move = moves{j} ;

        if isnumeric(this_move)
            i = this_move(1) ;

            isaccumulate = i ~= 1 || jj<burn_in ;
            
            % swap move if this_move has 2 indices
            if length(this_move) == 2 && isswap
                [swapX,swaps] = swap_randrect( X{i} , X{this_move(2)} , n_trials , ...
                                               flip_LL{i} , flip_LL{this_move(2)}) ;
                
                swapX.stats = stats{this_move(2)} ;

% figure(hswaps)
% imagesc( xor( swaps{1}.state>0 , swapX.state>0 ) )
% drawnow

                                           
                [ accumulated , swapX ] = ...
                    flip_MCMC( accumulated , swapX , accumulator , swaps , isaccumulate ) ;
                
                X{i} = swapX.X ; X{this_move(2)} = swapX.with ;
                stats{this_move(2)} = swapX.stats ;
                
            % regular MCMC move if this_move has one index
            elseif length(this_move) == 1
                i = this_move ;
                if jj>burn_in
                    jit_samples = jitter{i}(X{i},n_trials) ;
                else
                    jit_samples = jitter{i}(X{i},1) ;
                end                
                [ accumulated , X{i} ] = ...
                    flip_MCMC( accumulated , X{i} , accumulator , jit_samples , isaccumulate ) ;
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
            subplot(2,ceil(N_instances/2),i)
%             subplot(1,N_instances,i)
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
            titl = sprintf('X_%d   \\beta %.1f     excl %.2f',...
                            i,betas(i),exclusions(i)) ;
            if i > 1 && isswap
                s = stats{i} ;
                titl = sprintf('%s     acc %.2f',titl,s.accepted/s.N) ;
            end
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
cone_map.N_cones_prior  = func2str(N_cones_prior) ;
cone_map.betas          = betas ;
cone_map.n_cones        = n_cones ;
cone_map.burn_in        = burn_in ;
cone_map.N_iterations   = N_iterations ;
cone_map.moves          = moves ;
cone_map.accumulated    = reshape( accumulated , [M0 M1 N_colors] ) ;
cone_map.stats          = stats ;
cone_map.exclusions     = exclusions ;

figure
acc = cone_map.accumulated ./ max(cone_map.accumulated(:)) ;
imagesc(acc)
titl = sprintf('X_%d   \\beta %.2f    excl %.1f',i,betas(1),exclusions(1)) ;
if i == ceil(N_instances/4)
    title({sprintf('Iteration %d',jj) ; titl },'FontSize',22)
else
    title( titl , 'FontSize',22)
end


end