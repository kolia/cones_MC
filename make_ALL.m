% 
% % save('peach/greed_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters',...
% %      'greed_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')
% % save('peach/mcmc_peach_NROI2_ROI12_support3_beta02_delta01_SS4_5e05iters',...
% %      'mcmc_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')
% % save('peach/cast_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters',...
% %      'cast_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')

load('peach/greed_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')
load('peach/mcmc_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')
load('peach/cast_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')

type = 1 ;
main_CAST

this = pwd ; 
% try
make_plots( greed_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters,...
            mcmc_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters,...
            cast_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters, cone_map ) ; 
% catch e, cd(this) ; end

clear all

% try
    load('george/cast_george_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')
    load('george/mcmc_george_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')
    load('george/greed_george_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')
% catch
%     addpath('../agricola')
%     reap
% 
%     save('george/cast_george_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters',...
%          'cast_george_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')
% 
%     save('george/mcmc_george_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters',...
%          'mcmc_george_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')
% 
%     save('george/greed_george_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters',...
%         'cast_george_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')
% end

type = 0 ;
main_CAST

this = pwd ;
make_plots( greed_george_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters,...
            mcmc_george_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters,...
            cast_george_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters, cone_map ) ;
