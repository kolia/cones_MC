function cone_map = cones( cone_map )
%% cone_map = cones( LL )
%  Run MCMC to find cone locations.

% addpath('../libraries/swaps')

LL = cone_map.LL ;
N_cones_factor = cone_map.N_cones_factor ;

[M0,M1,N_colors] = size(LL) ;

%% plot and display info every ? MCMC iterations  (0 for never)
plot_every      = 20 ;
display_every   = 50 ;


%% PARAMETERS FOR MCMC
  TOTAL_trials  = 30 * M0 * M1 ;  % number of trials after burn-in = TOTAL_trials * n_trials ;
  burn_in       = 20 * M0 * M1 ; % number of burn-in trials

% maxcones      maximum number of cones allowed
  maxcones      = 150 ;
%   maxcones      = floor( 0.005 * M0 * M1 ) ;  
  
% these params shouldn't need tweaking unless the problem setup changes
%
% deltas        powers applied to likelihoods of instances
  deltas        = [0.4 0.35 0.3 0.2] ;
%   deltas        = make_deltas(0.1,0.4,4.4,20) ;

% betas         temperatures of independent instances run simultaneously
  betas         = ones(1,length(deltas)) ;

% exclusions    exclusion distance between cones of instance i
  exclusions    = ones(1,length(betas)) * 9.2  ;

% q             probability of trying to move an existing cone vs. placing
%               a new one.
  q             = 0.99 ;

  % moves         sequence of moves at each iteration, currently:
%               - a regular MC move for each instance


%% MCMC RUN
N_instances     = length(betas) ;
fprintf('\n\nSTARTING %d MCMC instances with',N_instances)
fprintf('\ndifferent inverse temperatures beta:\n')
fprintf('%.2f   ',betas)
fprintf('\ntemperature powers delta:\n')
fprintf('%.2f   ',deltas)
fprintf('\nmaximum number of cones: %d   ',maxcones)


% initializing variables
jitter          = cell( N_instances , 1 ) ;
X               = cell( N_instances , 1 ) ;
updater         = cell( N_instances , 1 ) ;
LLi             = cell( N_instances , 1 ) ;

for i=1:N_instances

	LLi{i}          = ((1+LL).^deltas(i)-1) * betas(i) ;
    N_factor_i      = ((1+N_cones_factor).^deltas(i)-1) * betas(i) ;
    
    jitter{i}       = @(X)move(X  , 2 , q , LLi{i} ) ;
  
    % initialize X{i}
    X{i}            = initialize_X( LLi{i} , N_factor_i , exclusions(i) , maxcones) ;
    
    updater{i}      = @(X,trial)update_X(X,trial,LLi{i}) ;
    
end


accumulator  = @(X)[(X.state(:)'==1) (X.state(:)'==2) (X.state(:)'==3)] ;
accumulated  = zeros( 1 , M0*M1*N_colors ) ;

moves = [num2cell(1:N_instances) num2cell([1:N_instances-1 ; 2:N_instances],1)] ;
% moves = num2cell(1:N_instances) ;  % no swaps for now

N_moves      = length(moves) ;
burn_in      = ceil( burn_in / N_moves ) ;
N_iterations = burn_in + ceil( TOTAL_trials / N_moves ) ;
n_cones      = zeros( N_iterations , 1 ) ;

stats = cell(N_moves,1) ;
for j=1:N_moves
    stats{j}.N500     = 0 ;
    stats{j}.accepted = 0 ;
end

if plot_every
scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)*0.7*0.5 1500*0.5 1200*0.5]) ;
end

% hswaps = figure ;

GG0 = LL - min(LL(:)) ;
GG0 = ( GG0 / max(GG0(:))) .^ 0.7 ; % !!! take power to enhance low values visually

% MAIN MCMC LOOP
% fprintf('\n\n      MCMC progress:|0%%              100%%|\n ')
fprintf('\n\nMCMC progress:')
t = cputime ;
tic
for jj=1:N_iterations
%     if ~mod(jj,floor(N_iterations/20)) , fprintf('*') , end

    isswap = jj>burn_in*0.6 ;

    for j=1:N_moves
        this_move = moves{j} ;

        if isnumeric(this_move)
            i = this_move(1) ;

            isaccumulate = i ~= 1 || jj<burn_in ;
            
            % swap move if this_move has 2 indices
            if length(this_move) == 2 && isswap
                [swapX,swaps] = swap_closure( X{i} , X{this_move(2)} , ...
                                               LLi{i} , LLi{this_move(2)}) ;

                swapX.stats = stats{this_move(2)} ;

% figure(hswaps)
% imagesc( xor( swaps{1}.state>0 , swapX.state>0 ) )
% drawnow
                                           
                [ accumulated , swapX ] = ...
                    flip_MCMC( accumulated , swapX , accumulator , swaps , ...
                    @(x,t)update_swap(x,t,LLi{i},LLi{this_move(2)}) , isaccumulate ) ;
                
                X{i} = swapX.X ; X{this_move(2)} = swapX.with ;
                stats{this_move(2)} = swapX.stats ;
                
            % regular MCMC move if this_move has one index
            elseif length(this_move) == 1
                i = this_move ;
                [ accumulated , X{i} ] = ...
                    flip_MCMC( accumulated ,     X{i} ,       accumulator , ...
                               jitter{i}(X{i}) , updater{i} , isaccumulate ) ;
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
        if isswap
            swapstring = [swapstring ' swap'] ;
        else
            swapstring = [swapstring '     '] ;
        end
        fprintf('\nIteration:%4d of %d \t %s\t  %4d cones \t%8.2f sec',...
                            jj,N_iterations,swapstring,n_cones(jj),toc)
        tic
    end
    
    if ~mod(jj,plot_every)
        figure(h)
        for i=1:N_instances
            NN = floor(sqrt(N_instances)) ;
            subplot(NN,ceil(N_instances/NN),i)
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
            titl = sprintf('\\delta %.2g',deltas(i)) ;
            if i > 1 && isswap
                s = stats{i} ;
                titl = sprintf('%s     acc %.2f',titl,s.accepted/s.N500) ;
            end
            if i == 2 %ceil(N_instances/8)
                title({sprintf('Iteration %d',jj) ; titl },'FontSize',16)
            else
                title( titl , 'FontSize',16)
            end
            %             set(get(gca,'Title'),'Visible','on')
        end
        drawnow
    end
    
end
fprintf('\ndone in %.1f sec\n\n',cputime - t) ;

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

cone_map.code           = file2str('cones.m') ;

cone_map.betas          = betas ;
cone_map.deltas         = deltas ;
cone_map.n_cones        = n_cones ;
cone_map.burn_in        = burn_in ;
cone_map.N_iterations   = N_iterations ;
cone_map.moves          = moves ;
cone_map.accumulated    = reshape( accumulated , [M0 M1 N_colors] ) ;
cone_map.stats          = stats ;
cone_map.exclusions     = exclusions ;

if plot_every
    figure
    acc = cone_map.accumulated ./ max(cone_map.accumulated(:)) ;
    imagesc(acc)
    titl = sprintf('X_%d   \\beta %.2g    excl %.1f',i,betas(1),exclusions(1)) ;
    if i == ceil(N_instances/4)
        title({sprintf('Iteration %d',jj) ; titl },'FontSize',16)
    else
        title( titl , 'FontSize',16)
    end
end

end