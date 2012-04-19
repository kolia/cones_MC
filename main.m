addpath(genpath(pwd))

%% LOAD DATA

% load george/stas       
load peach/peach_data 
% stas(i).spikes  : the list of spike times of cell i
% stas(i).spatial : the spatial component of the STA of cell i

% load george/cone_params   
load peach/cone_params    % contains struct cone_params with fields:
% stimulus_variance : variance of each stim pixel color chanel (usually 1)
% supersample       : the integer number of cone positions per pixel
%                     width/height  (usually 4) 
% colors            : 3x3 matrix of cone color sensitivities
% support_radius    : radius of cone receptive field filter  (usually 3.0)
% repulsion_radii


% restrict the stimulus to a region of interest  (if memory is a problem)
% stas = restrict_ROI( stas, [1 160], [1 160] ) ;


%% CONE_MAP contains parameters and preprocessed data structures

% calculate sparsity pattern (for MCMC) and initial LL map (for greedy)
cone_map = exact_LL_setup(stas,cone_params) ;

% How many MCMC/CAST iterations?
cone_map.N_iterations  = 1e4 ;
cone_map.max_time      = 4e5 ;  % in seconds

% % override defaults: plot, display, and save every ? iterations (0 = never)
% cone_map.plot_every  = 1000 ;
% cone_map.display_every = 20 ;
% cone_map.save_every    = 0  ;


%% run GREEDY, MCMC or CAST
% greed = greedy_cones(cone_map) ;  save('greed','greed') ;
% mcmc = MCMC(cone_map) ;           save('mcmc' ,'mcmc' ) ;
cast = CAST(cone_map) ;           save('cast' ,'cast' ) ;