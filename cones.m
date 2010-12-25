function [cone_map,results] = cones( cone_map , ID )
%% cone_map = cones( cone_map , ID )
%  Run MCMC to find cone locations. 

cone_map    = rmfield(cone_map,'ROI') ;

if nargin<2 , ID = 0 ; end
M0 = cone_map.M0 ;
M1 = cone_map.M1 ;

%% PARAMETERS FOR MCMC
  cone_map.N_iterations  = 10000 ; %6   * M0 * M1 ; % number of trials
  cone_map.start_swap    = 10 ; %1/2 * M0 * M1 ; % iteration w/ swapping starts

% betas         temperatures of independent instances run simultaneously
  cone_map.betas         = make_deltas(0.1,1,1,20) ;

% deltas        temperatures of independent instances run simultaneously
  cone_map.deltas        = make_deltas(0.75,1,2,length(cone_map.betas)) ;
  
% q             probability of trying to move an existing cone vs. placing
%               a new one.
%  cone_map.q             = 0.99 ;

% moves         sequence of moves at each iteration, currently:
%               - a regular MC move for each instance
% cone_map.moves = [num2cell(1:N_instances) ...
%                   num2cell(1:N_instances) ...
%                   num2cell([1:N_instances-1 ; 2:N_instances],1)] ;

[cone_map,results] = MCMC_parallel_tempering( cone_map , ID ) ;

end