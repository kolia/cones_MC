% PREPARE cone_map

warning off
type = 0 ;

if type==1
    type = 'peach' ;
    load peach/peach_data    % contains 'stas'
    load peach/cone_params
else
    type = 'george' ;
    load george/stas   % contains 'stas'
    load george/cone_params   % contains 'cone_params'
    cone_params.stimulus_variance = 1 ;
end

% ROIs = {[1 size(stas(1).spatial,1)] [1 size(stas(1).spatial,2)] [1 size(stas(1).spatial,1)/2] [1 size(stas(1).spatial,2)/2] [size(stas(1).spatial,1)/2+1 size(stas(1).spatial,1)] [size(stas(1).spatial,2)/2+1 size(stas(1).spatial,2)]} ;
% roi = 4 ;
% roj = 6 ;
ROIs = {[1 size(stas(1).spatial,1)] [1 size(stas(1).spatial,2)]} ;
roi = 1 ;
roj = 2 ;
stas = restrict_ROI( stas, ROIs{roi}, ROIs{roj} ) ;

cone_params.support_radius = 3 ;
% cone_params.supersample = 2 ;
cone_map = exact_LL_setup(stas,cone_params) ; % cone_map, aka PROB or data

imagesc(cone_map.NICE)

cone_map.N_iterations  = 5e5 ;
cone_map.max_time      = 2e5 ;
cone_map.betas  = make_deltas( 0.2, 1, 1, 20 ) ;
cone_map.deltas = make_deltas( 0.3, 1, 1, length(cone_map.betas) ) ;

cone_map.initX.rois   = [roi roj] ;
cone_map.initX.NROI   = numel(ROIs) ;
cone_map.initX.ROI    = {ROIs{roi} ; ROIs{roj}} ;
cone_map.initX.type   = type ;
cone_map.initX.supersample  = cone_params.supersample ;
cone_map.initX.support_radius = cone_params.support_radius ;
cone_map.initX.N_iterations   = cone_map.N_iterations ;
cone_map.initX.betas  = cone_map.betas  ;
cone_map.initX.deltas = cone_map.deltas ;

cone_map.plot_every    = 0 ;
cone_map.display_every = 20 ;

base_str = cone_map_string( cone_map ) ;

% % THEN RUN THIS to run on your own computer:
% greed = greedy_cones(cone_map) ;  save(['greed_' base_str],'greed')
% mcmc = MCMC(cone_map) ;           save(['mcmc_'  base_str],'mcmc' )
% cast = CAST(cone_map) ;           save(['cast_'  base_str],'cast' )

% % OR THIS to run 50 MCMC instances and 50 CAST on the hpc cluster:
% %            INSTALL AGRICOLA FIRST
N = 30 ;
ids = cell(1,N) ;
for i=1:length(ids) , ids{i} = {i} ; end

PBS.l.mem = '1500mb' ;
PBS.l.walltime = '70:00:00' ;
sow(['cast_' base_str],@(ID)CAST(cone_map,ID),ids,PBS) ;

PBS.l.mem = '1500mb' ;
PBS.l.walltime = '70:00:00' ;
sow(['mcmc_' base_str],@(ID)MCMC(cone_map,ID),ids,PBS) ;

sow(['greed_' base_str],@()greedy_cones(cone_map)) ;
