% PREPARE cone_map

% load peach/peach_data    % contains 'stas'
% load peach/cone_params

% load george/stas   % contains 'stas'
% load george/cone_params   % contains 'cone_params'

cone_params.fudge = 1 ;
cone_map = exact_LL_setup(stas,cone_params) ; % cone_map, aka PROB or data


cone_map.betas  = make_deltas( 0.2, 1, 1, 20 ) ;
cone_map.deltas = make_deltas( 0.3, 1, 1, length(cone_map.betas) ) ;

cone_map.plot_every    = 0  ;
cone_map.display_every = 20 ;

% % THEN RUN THIS to run on your own computer:
greed = greedy_cones(cone_map) ;
mcmc = MCMC(cone_map) ;
% cast = CAST(cone_map) ;

% OR THIS to run 50 MCMC instances and 50 CAST on the hpc cluster:
%            INSTALL AGRICOLA FIRST
% sow('greed',@()greedy_cones(cone_map)) ;
% N = 50 ;
% ids = cell(1,N) ;
% for i=1:length(ids) , ids{i} = {i} ; end
% sow('mcmc',@(ID)MCMC(cone_map,ID),ids) ;
% sow('cast',@(ID)CAST(cone_map,ID),ids) ;