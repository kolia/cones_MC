
% try
%     load('peach/greed_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')
%     load('peach/mcmc_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')
%     load('peach/cast_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')
% catch
%     addpath('../agricola')
%     reap
% 
%     save('peach/cast_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters',...
%          'cast_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')
% 
%     save('peach/mcmc_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters',...
%          'mcmc_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')
% 
%     save('peach/greed_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters',...
%         'cast_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters')
% end
% 
% type = 1 ;
% main_CAST
% 
% this = pwd ; 
% % try
% make_plots( greed_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters,...
%             mcmc_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters,...
%             cast_peach_NROI2_ROI12_support3_SS4_beta02_delta01_5e05iters, cone_map ) ; 
% % catch e, cd(this) ; end
% 
% clear all

% try
    load('george/cast_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters')
% %     for i=1:numel(cast_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters)
% %         cast_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters{i}.bestX = ...
% %             rmfield(cast_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters{i}.bestX,{'LL_history','cputime','N_cones_history','dX','excluded','sparse_STA_W_state','swap'}) ;
% %     end
    load('george/mcmc_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters')
%     for i=1:numel(mcmc_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters)
%         for j=1:numel(mcmc_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters{i}.bestX)
%             mcmc_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters{i}.bestX{j} = ...
%                 rmfield(mcmc_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters{i}.bestX{j},{'LL_history','cputime','N_cones_history','excluded','sparse_STA_W_state'}) ;
%         end
%     end
    load('george/greed_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters')
% catch
%     addpath('../agricola')
%     reap
% 
%     save('george/cast_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters',...
%          'cast_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters','-v7.3')
% 
%     save('george/mcmc_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters',...
%          'mcmc_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters')
% 
%     save('george/greed_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters',...
%         'greed_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters')
% end

type = 0 ;
main_CAST

this = pwd ;
make_plots( greed_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters,...
            mcmc_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters,...
            cast_george_NROI2_ROI12_support3_SS4_beta02_delta01_1e06iters, cone_map ) ;
