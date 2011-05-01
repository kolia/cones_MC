% load peach_data    % contains 'stas'
% load cone_params   % contains 'cone_params'
% 
% cone_map = exact_LL_setup(stas,cone_params) ; % cone_map, aka PROB or data

cone_map.betas  = make_deltas( 0.5,1,1,50) ;
cone_map.deltas = make_deltas(0.85,1,1,length(cone_map.betas)) ;

cone_map = CAST(cone_map) ;