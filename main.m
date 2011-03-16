function main( cone_map )

% if nargin>1
%     cone_map.SS = cone_map.cone_params.supersample ;
%     X = initialize_X(cone_map.M0,...
%                      cone_map.M1,...
%                      cone_map.N_colors,...
%                      cone_map.SS,...
%                      cone_map.cone_params.replusion_radii,1,1) ;
%     [x,y,c] = find( greedy.X.state ) ;
%     for i=1:numel(x)
%         X = change_cone( X , [x(i) y(i) c(i)] , cone_map ) ;
%         [dummy,X] = update_X({},{X X},2) ;
%     end
% end

% cone_map.betas          = make_deltas( 0.5,1,1,50) ;
% cone_map.deltas         = make_deltas(0.85,1,1,length(cone_map.betas)) ;
% cone_map.N_iterations   = Inf ;
% cone_map.max_time       = 3600 ;

N = 60 ;    % MULTIPLE OF 4, because of splits in INDS below !!!

ids = cell(1,N) ;
for i=1:length(ids) , ids{i} = {i} ; end

cone_map.max_time = 2000 ;
initID = sow( 'initMC' , @(ID)MCMC( cone_map , ID ), ids ) ;
pause( cone_map.max_time + 60 )


cone_map.max_time = 10000 ;
oneID  = sow( 'one', @(ID)MCMC_parallel_tempering( ...
        coalescer(cone_map, ['../' initID.id], 'bestX_[0-9]+.mat' , ...
        ID , @prep_cool , ID )), ids ) ;

pause( cone_map.max_time + 60 )
    
cone_map.max_time = 30000 ;

INDS =  { { 1 ; 1:floor(N/2) }  { 2 ; floor(N/2)+1:N }  { 3 ; 1:2:N } ...
        { 4 ; 2:2:N } {5 ; [(1:floor(N/4)) (floor(N/2)+1:floor(3*N/4))]} ...
        { 6 ; [(floor(N/4)+1:floor(N/2)) (floor(3*N/4):N)]} ...
        { 7 ; [1:4:N 2:4:N]} {8  ; [3:4:N 4:4:N]} ...
        { 9 ; [1:4:N 3:4:N]} {10 ; [2:4:N 4:4:N]}} ;

coolID = sow( 'cool', @(ID,is)MCMC_parallel_tempering( ...
        coalescer(cone_map, ['../' oneID.id], 'stats_[0-9]+.mat' , is , ...
        @prep_cool , ID )), INDS ) ;

% sow( 'PT', @(ID)MCMC_parallel_tempering( ...
%       coalescer(cone_map, ['../' coolID.id], 'stats_[0-9]+.mat' , ID , ...
%       @prep_PT , ID )), ids(1:length(INDS)) )

end



function cone_map = prep_cool( cone_map )

M0          = cone_map.M0 ;
M1          = cone_map.M1 ;
N_instances = length(cone_map.X) ;

%% PARAMETERS FOR MCMC
cone_map.N_iterations  = 1000 ; %1   * M0 * M1 ; % number of trials
cone_map.start_swap    = 0 ;

% betas         temperatures of independent instances run simultaneously
cone_map.betas         = make_deltas( 0.4, 1, 1, N_instances ) ;

% moves         sequence of moves at each iteration, currently:
%               - a regular MC move for each instance
% cone_map.moves = num2cell([ones(1,N_instances-1) ; 2:N_instances],1) ;  % al
cone_map.moves = [num2cell([1:N_instances-1 ; 2:N_instances],1)] ;    % line

cone_map.N_best = 1 ;

end


function cone_map = prep_PT( cone_map )

M0          = cone_map.M0 ;
M1          = cone_map.M1 ;
N_instances = 200 ;

%% PARAMETERS FOR MCMC
cone_map.N_iterations   = 20000 ;
cone_map.max_time       = 1000000 ;
cone_map.start_swap     = 10 ;

% betas         temperatures of independent instances run simultaneously
cone_map.betas         = make_deltas( 0.4, 1, 1, N_instances ) ;

% deltas        temperatures of independent instances run simultaneously
cone_map.deltas        = make_deltas(0.3,1,3,length(cone_map.betas)) ;
  
cone_map.N_best = 1 ;

end