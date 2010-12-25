function cone_map = exact_LL_setup( GC_stas , cone_params , cone_map )
%% cone_map = exact_LL_setup( GC_stas , cone_params , cone_map )
%  Expand data and parameters into variables used to calculate likelihoods.
%  Mainly, spatial supersampling by a factor of cone_map.supersample is
%  applied to the STAs, and the convolution of the STAs with the cone
%  receptive fields is stored in STA_W.

addpath(genpath(pwd))

if nargin < 3  ,   cone_map = struct ; end

% size of data
  [M0,M1,N_colors] = size(GC_stas(1).spatial) ;


%% PROBLEM SETUP: ganglion cell STAs

% supersample factor
SS  = cone_params.supersample ;
R   = ceil( cone_params.support_radius * SS ) ;

% set up Region of Interest
if ~isfield(cone_map,'ROI')
    x = repmat( 1/(2*SS):1/SS:M0-1/(2*SS) , 1 , M1*SS ) ;
    y = repmat( 1/(2*SS):1/SS:M1-1/(2*SS) , M0*SS , 1 ) ;
    cone_map.ROI = [x' y(:)] ;
    clear x y
end
NROI  = size(cone_map.ROI,1) ;

% Unpacking GC_stas into: STA, norms of STAs and N_spikes
N_GC = length(GC_stas) ;
STA_norm = zeros(N_GC,1) ;
N_spikes = zeros(N_GC,1) ;
STA      = zeros(M0*M1*N_colors,N_GC) ;
for i=1:N_GC
    N_spikes(i) = length(GC_stas(i).spikes) ;
    STA(:,i)    = GC_stas(i).spatial(:) ;
    STA_norm(i) = norm(STA(:,i)) ;
end

% cell_consts = N_spikes ./ exp(STA_norm/2) * cone_params.stimulus_variance ;
cell_consts = N_spikes * cone_params.stimulus_variance ;

% memoized(?) function returning gaussian mass in a box
gaus_in_box = gaus_in_a_box( cone_params.sigma , SS ) ;

prior_cov   = cone_params.stimulus_variance^2*N_GC/sum(STA_norm.^2) ;

cone_map.N_cones_term = sum( log( prior_cov) - log(cell_consts(:)+prior_cov) ) ;
cone_map.quad_factor  = N_spikes.^2 ./ (cell_consts+prior_cov) ;

% % stereotyped cone receptive field
% s = cone_params.sigma ;
% R = cone_params.support_radius * SS ;
% cone_RF = exp(-0.5 * ((-R:R)/(SS*s)).^2)' * exp(-0.5 * ((-R:R)/(SS*s)).^2) ;
% cone_RF = cone_RF / norm(cone_RF(:)) ; clear s



%% SETUP for Log-LIKELIHOOD calculations

try load('STA_W')  % calculating STA_W takes a while, check for STA_W.mat
catch
    
    STA_W = zeros(NROI*N_colors,N_GC) ;
    fprintf('Calculating STA_W...\n')
    for ii=1:NROI
        i   = cone_map.ROI(ii,1) ;
        j   = cone_map.ROI(ii,2) ;
        BW  = make_filter(M0,M1,i,j,cone_params.support_radius,gaus_in_box) ;
        
        for c=1:N_colors
            filter  = kron(cone_params.colors(c,:),BW(:)') ;
            STA_W((c-1)*NROI + ii , :) = filter * STA ;
        end
    end
    STA_W = STA_W ./ cone_params.stimulus_variance ; % repmat( (N_spikes ./ cell_consts)' ,NROI*N_colors,1) ;

save('STA_W','STA_W')
end

LL = (STA_W .^ 2) * diag(cell_consts) ;
LL = reshape( sum(LL , 2)/2 , [M0*SS M1*SS 3] ) ;
IC = inv(cone_params.colors) ;
QC = reshape( reshape(LL,[],3) * IC' , size(LL) ) ;
QC = QC - min(QC(:)) ;
NICE = QC ./ max(QC(:)) ;


coneConv = zeros( 2*R+SS , 2*R+SS , SS , SS ) ;

f = 1/(2*SS):1/SS:2*R/SS+1 ;
for xx=1:2*R+SS
    x = f(xx) ;
    for yy=1:2*R+SS
        y = f(yy) ;
        
        a = make_filter(4*R/SS+1,4*R/SS+1,x+SS,y+SS,cone_params.support_radius,gaus_in_box) ;
                
        for ss=1:SS
            s = (ss-0.5)/SS ;
            for tt=1:SS
                t = (tt-0.5)/SS ;
                
                b = make_filter(4*R/SS+1,4*R/SS+1,2*R/SS+s,2*R/SS+t,...
                    cone_params.support_radius,gaus_in_box) ;
                
                coneConv(xx,yy,ss,tt) = dot(a(:),b(:)) ;
            end
        end
    end
end


cone_map.R              = R ;
cone_map.coneConv       = coneConv ;
% cone_map.sumLconst      = sum(log(2*pi*cell_consts)) ;
cone_map.STA_W          = STA_W ;
cone_map.min_STA_W      = min(STA_W(:)) ;
cone_map.NICE           = NICE ;
cone_map.M0             = M0 ;
cone_map.M1             = M1 ;
cone_map.N_colors       = N_colors ;
cone_map.N_GC           = N_GC ;
cone_map.cell_consts    = cell_consts ;
cone_map.prior_cov      = prior_cov ;
cone_map.colorDot       = cone_params.colors * cone_params.colors' ;
cone_map.NROI           = NROI ;
cone_map.cone_params    = cone_params ;

end



function filter = make_filter(M0,M1,i,j,support,gaus_boxed)

filter = zeros(M0,M1) ;
ox  = max(1,floor(i-support)):min(M0,ceil(i+support)) ;
oy  = max(1,floor(j-support)):min(M1,ceil(j+support)) ;
x   = repmat(ox(:),numel(oy),1) ;
y   = reshape( repmat(oy,numel(ox),1) , [] , 1 ) ;
g   = gaus_boxed(i-x,j-y) ;
filter(x+(y-1)*M0) = g ;

end


function gf = gaus_in_a_box(sigma,SS)

% memo    = sparse([],[],[],1000,1000, ceil( (SS*3*sigma)^2 ) ) ;
gf      = @gib ;

    function out = gib(dx,dy)
        
%         out = zeros(size(in));  % preallocate output
%         [tf,loc] = ismember(in,x);  % find which in's already computed in x
%         ft = ~tf;  % ones to be computed
%         out(ft) = F(in(ft));  % get output values for ones not already in
%         % place new values in storage
%         x = [x in(ft(:).')];
%         y = [y reshape(out(ft),1,[])];
%         out(tf) = y(loc(tf));  % fill in the rest of the output values

        dx  = dx(:) ;
        dy  = dy(:) ;
        l   =  ones(size(dx)) ;
        O   = zeros(size(dx)) ;
        
        out = mvncdf([dx dy] + [l l],[],sigma.*[1 1]) ...
            - mvncdf([dx dy] + [O l],[],sigma.*[1 1]) ...
            - mvncdf([dx dy] + [l O],[],sigma.*[1 1]) ...
            + mvncdf([dx dy]        ,[],sigma.*[1 1]) ;
    end

end