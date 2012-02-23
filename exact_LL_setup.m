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

cone_params.support_radius = ceil(cone_params.support_radius) ;
cone_map.cone_params = cone_params ;

% set up Region of Interest
if ~isfield(cone_map,'ROI')
    x = repmat( 1/(2*SS):1/SS:M0-1/(2*SS) , 1 , M1*SS ) ;
    y = repmat( 1/(2*SS):1/SS:M1-1/(2*SS) , M0*SS , 1 ) ;
    ROI = [x' y(:)] ;
    clear x y
else
    ROI = cone_map.ROI ;
end
NROI  = size(ROI,1) ;

% Unpacking GC_stas into: STA, norms of STAs and N_spikes
N_GC = length(GC_stas) ;
STA_norm = zeros(N_GC,1) ;
N_spikes = zeros(N_GC,1) ;
STA      = zeros(M0,M1,N_colors,N_GC) ;
for i=1:N_GC
    N_spikes(i)  = length(GC_stas(i).spikes) ;
    STA(:,:,:,i) = reshape( GC_stas(i).spatial(:), M0, M1, N_colors ) ;
    STA_norm(i)  = norm(reshape(STA(:,:,:,i),1,[])) ;
end

% cell_consts = N_spikes ./ exp(STA_norm/2) * cone_params.stimulus_variance ;
cell_consts = N_spikes * cone_params.stimulus_variance ;

% memoized(?) function returning gaussian mass in a box
gaus_in_box_memo = gaus_in_a_box_memo( cone_params.sigma, SS, cone_params.support_radius ) ;

% prior_cov   = cone_params.stimulus_variance^2*N_GC/sum(STA_norm.^2) ;
% prior_cov   = cone_params.stimulus_variance^2*(N_GC-1)/sum(STA_norm.^2) ;
% prior_cov = cone_params.stimulus_variance^2*(N_GC-1)/sum(STA_norm.^2) ;

cone_map.prior_cov  = cone_params.stimulus_variance^2*...
                    (sum(N_spikes))/sum(N_spikes .* STA_norm.^2) ;
                
cone_map.cov_factor   = cell_consts+cone_map.prior_cov ;
cone_map.N_cones_term = sum( log( cone_map.prior_cov) - ...
                             log( cone_map.cov_factor(:)) ) ;
cone_map.quad_factor  = N_spikes.^2 ./ cone_map.cov_factor ;


cone_map.prior_covs    = (cone_params.stimulus_variance ./ STA_norm ).^2 ;
cone_map.cov_factors   = cell_consts+cone_map.prior_covs ;
cone_map.N_cones_terms = log( cone_map.prior_covs ) - log( cone_map.cov_factors) ;
cone_map.quad_factors  = N_spikes.^2 ./ cone_map.cov_factors ;


% sparse int matrix, with number of out-of-border adjacencies
cone_map.outofbounds = sparse([],[],[],M0*SS,M1*SS,2*(M0+M1)*SS) ;
cone_map.outofbounds(:,[1 M1*SS]) = 1 ;
cone_map.outofbounds([1 M0*SS],:) = cone_map.outofbounds([1 M0*SS],:) + 1 ;


%% SETUP for Log-LIKELIHOOD calculations

cone_map.M0             = M0 ;
cone_map.M1             = M1 ;
cone_map.N_colors       = N_colors ;
cone_map.N_GC           = N_GC ;
cone_map.cell_consts    = cell_consts ;
cone_map.NROI           = NROI ;
cone_map.N_spikes       = N_spikes ;
cone_map.SS             = SS ;

cone_map.naive_LL       = 0 ;

coneConv = zeros( 2*R+SS , 2*R+SS , SS , SS ) ;
WW = zeros(SS,SS) ;

f = 1/(2*SS):1/SS:2*R/SS+1 ;
for xx=1:2*R+SS
    x = f(xx) ;
    for yy=1:2*R+SS
        y = f(yy) ;
        
        a = make_filter_new(4*R/SS+1,4*R/SS+1,x+R/SS,y+R/SS,gaus_in_box_memo,...
                        cone_map.cone_params.support_radius) ;
        
        for ss=1:SS
            s = (ss-0.5)/SS ;
            for tt=1:SS
                t = (tt-0.5)/SS ;
                
                b = make_filter_new(4*R/SS+1,4*R/SS+1,2*R/SS+s,2*R/SS+t,gaus_in_box_memo,...
                                cone_map.cone_params.support_radius) ;
                
                coneConv(xx,yy,ss,tt) = dot(a(:),b(:)) ;
            end
        end
        if (xx<=SS) && (yy<=SS)
            WW(xx,yy) = dot(a(:),a(:)) ;
        end
    end
end

[cone_map.sparse_struct,cone_map.LL] = ...
    make_sparse_struct(cone_map,STA,WW,gaus_in_box_memo) ;

IC = inv(cone_params.colors) ;
QC = reshape( reshape(cone_map.LL+0.5*cone_map.N_cones_term,[],3) * IC', size(cone_map.LL) ) ;
% QC = reshape( reshape(cone_map.LL,[],3) * IC', size(cone_map.LL) ) ;
cone_map.NICE = plotable_evidence( QC ) ;

% imagesc( cone_map.NICE )


cone_map.R              = R ;
cone_map.gaus_boxed     = gaus_in_box_memo ;
cone_map.coneConv       = coneConv ;
cone_map.STA            = reshape( STA, [M0*M1*N_colors,N_GC] ) ;
cone_map.min_STA_W      = -0.2 ; %min(STA_W(:)) ;
cone_map.colorDot       = cone_params.colors * cone_params.colors' ;

% Some default values
cone_map.N_iterations   = 100000;
cone_map.q              = 0.5   ;
cone_map.plot_every     = 0     ;
cone_map.plot_skip      = 100   ;
cone_map.display_every  = 100   ;
cone_map.save_every     = 200   ;
cone_map.ID             = 0     ;
cone_map.max_time       = 200000;
cone_map.N_fast         = 1     ;

cone_map.initX = initialize_X( cone_map.M0, cone_map.M1, ...
                               cone_map.N_colors, cone_map.SS, ...
                               cone_map.cone_params.replusion_radii, ...
                               cone_map.naive_LL, 1, 1) ;

end

function [sparse_struct, LL] = make_sparse_struct(cone_map,STA,WW,gaus_boxed)

M0 = cone_map.M0 ;
M1 = cone_map.M1 ;
SS = cone_map.SS ;
N_colors = cone_map.N_colors ;
support = cone_map.cone_params.support_radius ;
colors  = cone_map.cone_params.colors ;

LL = zeros(M0*SS,M1*SS,N_colors) ;

supersamples = 1/(2*SS):1/SS:1 ;
gs = cell(SS) ;
sparse_struct = struct('x', cell(cone_map.N_GC,1), ...
              'y', cell(cone_map.N_GC,1), ...
              'c', cell(cone_map.N_GC,1)) ;

for ii=1:SS
    for jj=1:SS
        i = supersamples(ii) ;
        j = supersamples(jj) ;
        g = reshape( gaus_boxed(i,j), [2*support+1 2*support+1]) ;
        gs{ii,jj} = g(end:-1:1,end:-1:1) ;
    end
end

fprintf('making sparse struct and LL for GC number')
for gc=1:cone_map.N_GC
    gcLL = zeros(M0*SS,M1*SS,N_colors) ;
    for ii=1:SS
        for jj=1:SS    
            CC = zeros(M0*M1,N_colors) ;
            for color=1:N_colors
                CCC = conv2( STA(:,:,color,gc), gs{ii,jj} ) ;
                CCC = CCC(support+1:M0+support,support+1:M1+support) ;
                CC(:,color) = CCC(:) ;
            end
            C = 0.5 * cone_map.quad_factor(gc) * (CC * colors').^2 / WW(ii,jj) ;
            gcLL( ii:SS:M0*SS, jj:SS:M1*SS, :) = ...
                gcLL( ii:SS:M0*SS, jj:SS:M1*SS, :) + reshape(C,[M0 M1 3]) ;
        end
    end
    [sparse_struct(gc).x,sparse_struct(gc).y,sparse_struct(gc).c] = ...
            find( gcLL+cone_map.N_cones_terms(gc)>0 ) ;
    LL = LL + gcLL ;
    fprintf(' %d',gc)
end
end


function filter = make_filter_new(M0,M1,i,j,gaus_boxed, support)
filter = zeros(M0,M1) ;
[g,index] = filter_index( i, j, M0, M1, gaus_boxed, support) ;
filter(index) = g ;
end