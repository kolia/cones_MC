% add subfolders to path
addpath(genpath(pwd))


%% LOAD DATA

load george/stas       
% load peach/peach_data 
% stas(i).spikes  : the list of spike times of cell i
% stas(i).spatial : the spatial component of the STA of cell i

load george/cone_params   
% load peach/cone_params   
% cone_params.stimulus_variance : variance of each stim pixel color chanel
%                                 (usually 1)
%
% cone_params.supersample       : the integer number of cone positions per
%                                 pixel width/height
%                                 (usually 4)
%
% cone_params.colors            : 3x3 matrix of cone color sensitivities
%
% cone_params.support_radius    : radius of cone receptive field filter
%                                 (usually 3.0)
% 
% cone_params.repulsion_radii


% restrict the stimulus to a region of interest  (if memory is a problem)
% stas = restrict_ROI( stas, [1 160], [1 160] ) ;



%% CONE_MAP contains parameters and preprocessed data structures

% calculate sparsity pattern (for MCMC) and initial LL map (for greedy)
cone_map = exact_LL_setup(stas,cone_params) ;

% How many MCMC/CAST iterations?
cone_map.N_iterations  = 1e4 ;
cone_map.max_time      = 4e5 ;  % in seconds

% Define the progression of inverse temperatures for CAST
cone_map.min_delta = 0.1 ;
cone_map.min_beta  = 0.2 ;
cone_map.betas  = make_deltas( cone_map.min_beta, 1, 1, 20 ) ;
cone_map.deltas = make_deltas( cone_map.min_delta, 1, 1, length(cone_map.betas) ) ;

% plot, display, and save every N iterations (0 = never)
cone_map.plot_every    = 0   ;
cone_map.display_every = 20  ;
cone_map.save_every    = 0   ;



%% GREEDY
% greed = greedy_cones(cone_map) ;  save('greed','greed') ;



%% MCMC and CAST initialization

% track contacts during 'hot' greedy method
cone_map.track_contacts = true ;

% run 'hot' greedy method
greed_hot = greedy_cones(cone_map, 'hot') ;

% initialize MCMC and CAST with 'hot' greedy configuration
cone_map.initX = greed_hot.X ;

% plots are nice for MCMC and CAST
cone_map.plot_every    = 100   ;


%% MCMC or CAST

% % run MCMC
% mcmc = MCMC(cone_map) ;           save('mcmc' ,'mcmc' ) ;

% run CAST
cast = CAST(cone_map) ;           save('cast' ,'cast' ) ;