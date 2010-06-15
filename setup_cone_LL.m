function cone_map = setup_cone_LL( GC_stas , cone_params , cone_map )
%% cone_map = setup_cone_LL( GC_stas , cone_params , cone_map )
%  Expand data and parameters into variables used to calculate likelihoods.
%  Mainly, spatial supersampling by a factor of cone_map.supersample is
%  applied to the STAs, and the convolution of the STAs with the cone
%  receptive fields is stored in STA_W.

addpath(genpath(pwd))

if nargin < 3  ,   cone_map = struct ; end

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
R = cone_params.support_radius * SS ;
cone_RF = exp(-0.5 * ((-R:R)/(SS*s)).^2)' * exp(-0.5 * ((-R:R)/(SS*s)).^2) ;
cone_RF = cone_RF / norm(cone_RF(:)) ; clear s

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
%     BW      = imfilter(BW,cone_RF) ;
    
    for c=1:N_colors
        filter  = kron(cone_params.colors(c,:),BW(:)') ;
        STA_W((c-1)*NROI/N_colors + ii , :) = filter * STA ;
    end
end
STA_W = STA_W .* repmat( (N_spikes ./ cell_consts)' ,NROI,1) ;

save('STA_W','STA_W')
end

LL = (STA_W .^ 2) * diag(cell_consts) ;
LL = reshape( sum(LL , 2)/2 , [M0*SS M1*SS 3] ) ;

coneConv = zeros( 2*R+SS , 2*R+SS , SS , SS ) ;
for x=1:2*R+SS
    for y=1:2*R+SS
        a = placeRF(4*R+SS,cone_RF,R+x,R+y,SS) ;
        for s=1:SS
            for t=1:SS
                b = placeRF(4*R+SS,cone_RF,2*R+s,2*R+t,SS) ;
                v = dot(a(:),b(:)) ;
                coneConv(x,y,s,t) = v ;
            end
        end
    end
end

cone_map.R              = R ;
cone_map.coneConv       = coneConv ;

cone_map.sumLconst      = length(cell_consts) * log(2*pi) + ...
                          sum(log(cell_consts)) ;


cone_map.STA_W          = STA_W ;
cone_map.N_cones_factor = sum(log(2*pi*cell_consts)) ;
cone_map.LL             = LL ;
cone_map.M0             = M0 ;
cone_map.M1             = M1 ;
cone_map.M2             = M2 ;
cone_map.N_colors       = N_colors ;
cone_map.N_GC           = N_GC ;
cone_map.max_overlap    = coneConv(ceil(size(coneConv,1)/2),ceil(size(coneConv,2)/2)) ;
cone_map.cell_consts    = cell_consts ;
cone_map.colorDot       = cone_params.colors * cone_params.colors' ;
cone_map.ROIlogic       = ROIlogic ;
cone_map.ROI            = ROI ;
cone_map.N              = N ;
cone_map.NROI           = NROI ;
cone_map.SS             = SS ;

end


function a = placeRF(N,RF,x,y,SS)

R = floor(size(RF,1)/2) ;
a = zeros(N) ;
a(x-R:x+R,y-R:y+R) = RF ;
a = sum( reshape(a,SS,[]) ) ;
a = reshape(a,N/SS,[])' ;
a = sum( reshape(a,SS,[]) ) ;
a = reshape(a',N/SS,N/SS) ;
a = a/SS^2 ;
end