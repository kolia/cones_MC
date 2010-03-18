function [STA_W,cone_map] = setup_cone_LL( GC_stas , cone_params , cone_map )
%% [STA_W,cone_map] = setup_cone_LL( GC_stas , cone_params , cone_map )
%  Expand data and parameters into variables used to calculate likelihoods.
%  Mainly, spatial supersampling by a factor of cone_map.supersample is
%  applied to the STAs, and the convolution of the STAs with the cone
%  receptive fields is stored in STA_W.

addpath(genpath(pwd))

% size of data
  [M0,M1,N_colors] = size(GC_stas(1).spatial) ;
  M2 = M0*M1 ;


%% PROBLEM SETUP: ganglion cell STAs

% supersample factor
SS  = cone_params.supersample ;

% set up Region of Interest
if isfield(cone_map,'ROI')
    ROIlogic = zeros(M0*SS,M1*SS,N_colors) ;
    for c=1:N_colors
        ROIlogic(:,:,c) = kron(cone_map.ROI(:,:,c),ones(SS,SS)) ;
    end
else
    ROIlogic = ones(M0*SS,M1*SS,N_colors) ;
end
ROI = find(ROIlogic) ;
[I,J,dummy]   = ind2sub([M0*SS M1*SS N_colors],ROI) ;
sizeROI = [max(I)-min(I) max(J)-min(J)] + 1 ; 


% stereotyped cone receptive field
s = cone_params.sigma ;
r = cone_params.support_radius ;
cone_RF = exp(-0.5 * ((-SS*r:SS*r)/(SS*s)).^2)' * exp(-0.5 * ((-SS*r:SS*r)/(SS*s)).^2) ;
cone_RF = cone_RF / sum(cone_RF(:)) ; clear r s

N       = M2*(SS^2)*N_colors ;
NROI    = sum(ROI(:)>0) ;


% Unpacking GC_stas into: STA, norms of STAs and N_spikes
N_GC = length(GC_stas) ;
STA_norm = zeros(N_GC,1) ;
N_spikes = zeros(N_GC,1) ;
STA      = zeros(N,N_GC) ;
for i=1:N_GC
    N_spikes(i) = length(GC_stas(i).spikes) ;
    temp = cell(N_colors,1) ;
    for c=1:N_colors
        temp{c} = kron(GC_stas(i).spatial(:,:,c),ones(SS,SS)) ;
        temp{c} = temp{c}(:) ;
    end
    STA(:,i)    = cell2mat(temp) ;
    STA_norm(i) = norm(STA(:,i)) ;
end

cell_consts = N_spikes ./ exp(STA_norm/2) * cone_params.stimulus_variance ;

%% SETUP for Log-LIKELIHOOD calculations

try load('STA_W')  % calculating STA_W takes a while, check for STA_W.mat
catch
    
STA_W = zeros(NROI,N_GC) ;
fprintf('Calculating STA_W...\n')
for ii=1:NROI/N_colors
    i = I(ii) ;
    j = J(ii) ;
    BW      = zeros(M0*SS,M1*SS) ;
    BW(i,j) = 1 ;
    BW      = imfilter(BW,cone_RF) ;
    
    for c=1:N_colors
        filter  = kron(cone_params.colors(c,:),BW(:)') ;
        
        STA_W((c-1)*NROI/N_colors + ii , :) = filter * STA ;
    end
end
STA_W = STA_W .* repmat( (N_spikes ./ cell_consts)' ,NROI,1) ;

save('STA_W','STA_W')
end

coneConv    = conv2(cone_RF,cone_RF) ;

cone_map.M0             = M0 ;
cone_map.M1             = M1 ;
cone_map.M2             = M2 ;
cone_map.N_colors       = N_colors ;
cone_map.N_GC           = N_GC ;
cone_map.max_overlap    = coneConv(ceil(size(coneConv,1)/2),ceil(size(coneConv,2)/2)) ;
cone_map.cell_consts    = cell_consts ;
cone_map.stas           = GC_stas ;
cone_map.cone_RF        = cone_RF ;
cone_map.cone_params    = cone_params ;
cone_map.coneConv       = coneConv ;
cone_map.colorDot       = cone_params.colors * cone_params.colors' ;
cone_map.ROIlogic       = ROIlogic ;
cone_map.ROI            = ROI ;
cone_map.sizeROI        = sizeROI ;
cone_map.N              = N ;
cone_map.NROI           = NROI ;
cone_map.SS             = SS ;
