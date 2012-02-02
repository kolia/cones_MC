% PREPARE cone_map

type = 1 ;

if type==0
    type = 'peach' ;
    load peach/peach_data    % contains 'stas'
    load peach/cone_params
else
    type = 'george' ;
    load george/stas   % contains 'stas'
    load george/cone_params   % contains 'cone_params'
end

ROIs = {[1 42] [38 82] [78 122] [118 160]} ;

roi = 3 ;
roj = 3 ;
stas = restrict_ROI( stas, ROIs{roi}, ROIs{roj} ) ;

cone_params.fudge = 2 ;
cone_params.support_radius = 3 ;
cone_map = exact_LL_setup(stas,cone_params) ; % cone_map, aka PROB or data

cone_map.betas  = make_deltas( 0.2, 1, 1, 20 ) ;
cone_map.deltas = make_deltas( 0.3, 1, 1, length(cone_map.betas) ) ;

cone_map.plot_every    = 0  ;
cone_map.display_every = 20 ;

base_str  = sprintf('_%s_ROI-%d-%d_support-%d_fudge-%d',type,roi,roj,...
                     cone_params.support_radius,cone_params.fudge) ;

% % THEN RUN THIS to run on your own computer:
% greed = greedy_cones(cone_map) ;  save(['greed' base_str],'greed')
% mcmc = MCMC(cone_map) ;           save(['mcmc'  base_str],'mcmc' )
cast = CAST(cone_map) ;           save(['cast'  base_str],'cast' )

% OR THIS to run 50 MCMC instances and 50 CAST on the hpc cluster:
%            INSTALL AGRICOLA FIRST
% sow(['greed' base_str],@()greedy_cones(cone_map)) ;
% N = 30 ;
ids = cell(1,N) ;
% for i=1:length(ids) , ids{i} = {i} ; end
% 
% PBS.l.mem = '2000mb' ;
% PBS.l.walltime = '48:00:00' ;
% sow(['mcmc' base_str],@(ID)MCMC(cone_map,ID),ids,PBS) ;
% 
% PBS.l.mem = '3000mb' ;
% PBS.l.walltime = '72:00:00' ;
% sow(['cast' base_str],@(ID)CAST(cone_map,ID),ids,PBS) ;