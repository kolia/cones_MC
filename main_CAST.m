type = 1 ;

warning off

if type==1
    type = 'peach' ;
    load peach/peach_data     % contains 'stas'
    load peach/cone_params
else
    type = 'george' ;
    load george/stas          % contains 'stas'
    load george/cone_params   % contains 'cone_params'
    cone_params.stimulus_variance = 1 ;
end

cone_map = exact_LL_setup(stas,cone_params) ; % cone_map, aka PROB or data

cone_map.N_iterations  = 1e3 ;
cone_map.max_time      = 4e5 ;

cone_map.initX.type   = type ;
cone_map.type         = type ;
cone_map.initX.N_iterations   = cone_map.N_iterations ;

cone_map.plot_every    = 0 ;
base_str = cone_map_string( cone_map ) ;

% THEN RUN THIS to run on your own computer:
greed = greedy_cones(cone_map) ;  save(['greed_' base_str],'greed') ;
mcmc = MCMC(cone_map) ;           save(['mcmc_'  base_str],'mcmc' )
cast = CAST(cone_map) ;           save(['cast_'  base_str],'cast' )

% OR THIS to run 50 MCMC instances and 50 CAST on the hpc cluster:
%            INSTALL AGRICOLA FIRST
% N = 30 ;
% ids = cell(1,N) ;
% for i=1:length(ids) , ids{i} = {i} ; end
% 
% PBS.l.mem = '1500mb' ;
% PBS.l.walltime = '70:00:00' ;
% sow(['cast_' base_str],@(ID)CAST(cone_map,ID),ids,PBS) ;
% 
% PBS.l.mem = '1500mb' ;
% PBS.l.walltime = '70:00:00' ;
% sow(['mcmc_' base_str],@(ID)MCMC(cone_map,ID),ids,PBS) ;
% 
% sow(['greed_' base_str],@()greedy_cones(cone_map)) ;
