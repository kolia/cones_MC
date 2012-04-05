function cone_map = make_cone_map( type )

% PREPARE cone_map

warning off

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

ROIs = {[1 size(stas(1).spatial,1)] [1 size(stas(1).spatial,2)]} ;
roi = 1 ;
roj = 2 ;
stas = restrict_ROI( stas, ROIs{roi}, ROIs{roj} ) ;

cone_params.support_radius = 3 ;
% cone_params.supersample = 2 ;
cone_map = exact_LL_setup(stas,cone_params) ; % cone_map, aka PROB or data

imagesc(cone_map.NICE)

cone_map.N_iterations  = 5e5 ;
cone_map.max_time      = 4e5 ;
cone_map.profile_every = 0 ;
cone_map.min_delta = 0.1 ;
cone_map.min_beta  = 0.2 ;
cone_map.betas  = make_deltas( cone_map.min_beta, 1, 1, 20 ) ;
cone_map.deltas = make_deltas( cone_map.min_delta, 1, 1, length(cone_map.betas) ) ;

cone_map.initX.rois   = [roi roj] ;
cone_map.rois         = [roi roj] ;
cone_map.initX.NROIs  = numel(ROIs) ;
cone_map.NROIs        = numel(ROIs) ;
cone_map.initX.ROI    = {ROIs{roi} ; ROIs{roj}} ;
cone_map.initX.type   = type ;
cone_map.type         = type ;
cone_map.initX.supersample  = cone_params.supersample ;
cone_map.initX.support_radius = cone_params.support_radius ;
cone_map.initX.N_iterations   = cone_map.N_iterations ;
cone_map.initX.betas  = cone_map.betas  ;
cone_map.initX.deltas = cone_map.deltas ;

cone_map.plot_every    = 0 ;
cone_map.display_every = 20 ;