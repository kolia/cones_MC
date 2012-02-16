function base_str = cone_map_string(cone_map)

log_iters = round( log(cone_map.N_iterations)/log(10) ) ;

base_str  = sprintf('%s_NROI%d_ROI%d%d_support%d_fudge%d_SS%d_%de%diters',...
                     cone_map.initX.type,cone_map.initX.NROI,...
                     cone_map.initX.rois(1),cone_map.initX.rois(2),...
                     cone_map.cone_params.support_radius,...
                     cone_map.cone_params.fudge,...
                     cone_map.cone_params.supersample,...
                     round(cone_map.N_iterations/10^log_iters),log_iters) ;

end