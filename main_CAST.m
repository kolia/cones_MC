% PREPARE cone_map

load peach_data    % contains 'stas'
load cone_params   % contains 'cone_params'

cone_map = exact_LL_setup(stas,cone_params) ; % cone_map, aka PROB or data

cone_map.betas  = make_deltas( 0.2,1,1,20) ;
cone_map.deltas = make_deltas(0.3,1,1,length(cone_map.betas)) ;

cone_map.plot_every    = 20 ;
cone_map.display_every = 20 ;

cone_map.N_fast        = 3 ;

% THEN RUN THIS:
% cone_map = CAST(cone_map) ;