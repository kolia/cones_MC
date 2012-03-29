function cone_map = remake_cone_map( X )

if strcmp(X.type,'peach')
    type = 'peach' ;
    load peach/peach_data    % contains 'stas'
    load peach/cone_params
else
    type = 'george' ;
    load george/stas   % contains 'stas'
    load george/cone_params   % contains 'cone_params'
end

stas = restrict_ROI( stas, X.ROI{1}, X.ROI{2} ) ;

cone_params.support_radius = X.support_radius ;
cone_map = exact_LL_setup(stas,cone_params) ; % cone_map, aka PROB or data

cone_map.initX.rois   = X.rois  ;
cone_map.rois         = X.rois  ;
cone_map.initX.NROI   = X.NROI  ;
cone_map.initX.NROIs  = numel(X.ROI) ;
cone_map.NROIs        = numel(X.ROI) ;
cone_map.initX.ROI    = X.ROI   ;
cone_map.initX.type   = X.type  ;
cone_map.type         = X.type  ;
cone_map.initX.supersample = X.supersample ;
cone_map.initX.support_radius = X.support_radius ;
cone_map.min_beta     = min(X.betas) ;
cone_map.min_delta    = min(X.deltas) ;
cone_map.initX.betas  = X.betas  ;
cone_map.initX.deltas = X.deltas ;
cone_map.initX.N_iterations  = X.N_iterations ;

end