
% PREPARE cone_map

load george/stas       
% stas(i).spikes  : the list of spike times of cell i
% stas(i).spatial : the spatial component of the STA of cell i

load george/cone_params   
% cone_params.stimulus_variance : variance of each stim pixel color chanel
%                                 (default 1)
%
% cone_params.supersample       : the integer number of cone positions per
%                                 pixel width/height
%                                 (default 4)
%
% cone_params.colors            : 3x3 matrix of cone color sensitivities
%
% cone_params.support_radius    : radius of cone receptive field filter
%                                 (default 3.0)
% 
% cone_params.repulsion_radii   : 

stas = restrict_ROI( stas, [1 160], [1 160] ) ;

cone_params.support_radius = 3 ;
cone_map = exact_LL_setup(stas,cone_params) ;

cone_map.N_iterations  = 1e6 ;
cone_map.max_time      = 4e5 ;  % in seconds
cone_map.min_delta = 0.1 ;
cone_map.min_beta  = 0.2 ;
cone_map.betas  = make_deltas( cone_map.min_beta, 1, 1, 20 ) ;
cone_map.deltas = make_deltas( cone_map.min_delta, 1, 1, length(cone_map.betas) ) ;

% cone_map.initX.supersample  = cone_params.supersample ;
% cone_map.initX.support_radius = cone_params.support_radius ;
% cone_map.initX.N_iterations   = cone_map.N_iterations ;
% cone_map.initX.betas  = cone_map.betas  ;
% cone_map.initX.deltas = cone_map.deltas ;

cone_map.plot_every    =  0 ;   % zero means do not plot
cone_map.display_every = 20 ;   % display info every 20 iterations

greed_fast = greedy_cones(cone_map, 'fast') ;
cone_map.initX = greed_fast.X ;

greed = greedy_cones(cone_map) ;  save('greed','greed') ;
mcmc = MCMC(cone_map) ;           save('mcmc' ,'mcmc' ) ;
cast = CAST(cone_map) ;           save('cast' ,'cast' ) ;