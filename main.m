cone_map.datafolder = 'peach' ;
addpath(genpath(pwd))

%% LOAD DATA

load(sprintf('%s/stas', cone_map.datafolder))  % stas struct contains fields:
% stas(i).spikes    : the list of spike times of cell i
% stas(i).spatial   : the spatial component of the STA of cell i

load(sprintf('%s/cone_params', cone_map.datafolder)) % cone_params struct fields:
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
cone_map = exact_LL_setup(stas,cone_params,cone_map) ;

% How many MCMC/CAST iterations?
cone_map.N_iterations  = 5e5 ;
cone_map.max_time      = 4e5 ;  % in seconds

% % override defaults: plot, display, and save every ? iterations (0 = never)
% cone_map.plot_every  = 1000 ;
% cone_map.display_every = 20 ;
% cone_map.save_every    = 0  ;

% string identifying this cone_map problem; useful for naming saved results
cone_map.description = cone_map_string( cone_map ) ;

%% run GREEDY, MCMC or CAST locally
% result = GREEDY(cone_map) ;  save(['greed_' cone_map.description],'result') ;
result = MCMC(cone_map) ;           save(['mcmc_'  cone_map.description], 'result') ;
% result = CAST(cone_map) ;           save(['cast_'  cone_map.description], 'result') ;

%% some plots
make_plots( result )

%% %% or RUN 1 greedy instance, N MCMC instances and N CAST on the hpc cluster:
% %  start by cloning AGRICOLA from https://github.com/kolia/agricola
% addpath('../agricola')
% N = 30 ;  ids = cell(1,N) ; for i=1:length(ids) , ids{i} = {i} ; end
% cone_map.save_disk_space = true ;  % strip results of large fields
% 
% PBS.l.mem = '1500mb' ;
% PBS.l.walltime = '70:00:00' ;
% sow(['cast_'  base_str],@(ID)CAST(cone_map,ID),ids,PBS) ;
% sow(['mcmc_'  base_str],@(ID)MCMC(cone_map,ID),ids,PBS) ;
% sow(['greed_' base_str],@(  )GREEDY(cone_map)) ;