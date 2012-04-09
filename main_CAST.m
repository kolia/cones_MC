% function main_CAST( type )
type = 0 ;

cone_map = make_cone_map( type ) ;

base_str = cone_map_string( cone_map ) ;

% THEN RUN THIS to run on your own computer:
% try   
%     load(['greed_' base_str])
% catch e
%     greed = greedy_cones(cone_map) ;  save(['greed_' base_str],'greed') ;
% end

greed_fast = greedy_cones(cone_map, 'fast') ;
cone_map.initX = greed_fast.X ;

% mcmc = MCMC(cone_map) ;           save(['mcmc_'  base_str],'mcmc' )
% cast = CAST(cone_map) ;           save(['cast_'  base_str],'cast' )

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
