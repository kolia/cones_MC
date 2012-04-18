function cone_map = exact_LL_setup( GC_stas , cone_params , cone_map )
%% cone_map = exact_LL_setup( GC_stas , cone_params , cone_map )
%  Expand data and parameters into variables used to calculate likelihoods.
%  Mainly, spatial supersampling by a factor of cone_map.supersample is
%  applied to the STAs, and the convolution of the STAs with the cone
%  receptive fields is used to calculate the sparsity structure of
%  connections between cones and GCs, as well as a map of LL increases.

if nargin < 3  ,   cone_map = struct ; end

% size of data
[M0,M1,cone_map.N_colors] = size(GC_stas(1).spatial) ;


%% unpack data from GC_stas

% supersample factor
SS  = cone_params.supersample ;

% radius used for convolutions with cone RFs
R   = ceil( cone_params.support_radius * SS ) ;
cone_params.support_radius = ceil(cone_params.support_radius) ;

% copy cone_params
cone_map.cone_params = cone_params ;

% set up Region of Interest if it doesn't already exist
% these are fractional (x,y) coordinates, due to supersampling
if ~isfield(cone_map,'ROI')
    x = repmat( 1/(2*SS):1/SS:M0-1/(2*SS) , 1 , M1*SS ) ;
    y = repmat( 1/(2*SS):1/SS:M1-1/(2*SS) , M0*SS , 1 ) ;
    ROI = [x' y(:)] ;
    clear x y
else
    ROI = cone_map.ROI ;
end

% total number of possible cone positions
cone_map.NROI  = size(ROI,1) ;

% unpack GC_stas into: supersampled STAs, norms of STAs and N_spikes
N_GC = length(GC_stas) ;
STA_norm = zeros(N_GC,1) ;
cone_map.N_spikes = zeros(N_GC,1) ;
STA      = zeros(cone_map.N_colors,M0,M1,N_GC) ;
for i=1:N_GC
    cone_map.N_spikes(i)  = length(GC_stas(i).spikes) ;
    STA(:,:,:,i) = reshape(permute(GC_stas(i).spatial,[3 1 2]), cone_map.N_colors, M0, M1 ) ;
    STA_norm(i)  = norm(reshape(STA(:,:,:,i),1,[])) ;
end

%% calculate some constants etc in advance, to speed up actual calculations

% memoized function returning gaussian mass in a square pixel box
cone_map.gaus_boxed = gaus_in_a_box_memo( cone_params.sigma, SS, cone_params.support_radius ) ;

% constants used to calculate the log-likelihood
cone_map.cell_consts   = cone_map.N_spikes * cone_params.stimulus_variance ;
cone_map.prior_covs    = (cone_params.stimulus_variance ./ STA_norm ).^2 ;
cone_map.cov_factors   = cone_map.cell_consts+cone_map.prior_covs ;
cone_map.N_cones_terms = log( cone_map.prior_covs ) - log( cone_map.cov_factors) ;
cone_map.quad_factors  = cone_map.N_spikes.^2 ./ cone_map.cov_factors ;

% for fast lookup of legal shift moves in move.m
cone_map.outofbounds = sparse([],[],[],M0*SS,M1*SS,2*(M0+M1)*SS) ;
cone_map.outofbounds(:,[1 M1*SS]) = 1 ;
cone_map.outofbounds([1 M0*SS],:) = cone_map.outofbounds([1 M0*SS],:) + 1 ;

% copy params
cone_map.M0             = M0 ;
cone_map.M1             = M1 ;
cone_map.N_GC           = N_GC ;
cone_map.SS             = SS ;

%% make lookup table for dot products of cone RFs in all possible positions
cone_map.coneConv = zeros( 2*R+SS , 2*R+SS , SS , SS ) ;
WW = zeros(SS,SS) ;

f = 1/(2*SS):1/SS:2*R/SS+1 ;
for xx=1:2*R+SS
    x = f(xx) ;
    for yy=1:2*R+SS
        y = f(yy) ;
        
        a = make_filter_new(4*R/SS+1,4*R/SS+1,x+R/SS,y+R/SS,cone_map.gaus_boxed,...
                        cone_map.cone_params.support_radius) ;
        
        for ss=1:SS
            s = (ss-0.5)/SS ;
            for tt=1:SS
                t = (tt-0.5)/SS ;
                
                b = make_filter_new(4*R/SS+1,4*R/SS+1,2*R/SS+s,2*R/SS+t,...
                       cone_map.gaus_boxed,cone_map.cone_params.support_radius) ;
                
                cone_map.coneConv(xx,yy,ss,tt) = dot(a(:),b(:)) ;
            end
        end
        if (xx<=SS) && (yy<=SS)
            WW(xx,yy) = dot(a(:),a(:)) ;
        end
    end
end

% calculate sparsity structure and map of LL increases
[cone_map.sparse_struct,cone_map.LL] = ...
    make_sparse_struct(cone_map,STA,WW,gaus_in_box_memo) ;

% 
IC = inv(cone_params.colors) ;
QC = reshape( reshape(cone_map.LL,[],3) * IC', size(cone_map.LL) ) ;
cone_map.NICE = plotable_evidence( QC ) ;

% imagesc( cone_map.NICE )


cone_map.R              = R ;
cone_map.STA            = single( STA ) ;
cone_map.min_STA_W      = -0.2 ; %min(STA_W(:)) ;
cone_map.colorDot       = cone_params.colors * cone_params.colors' ;

% Some default values
cone_map.N_iterations   = 100000;
cone_map.q              = 0.5   ;
cone_map.plot_every     = 0     ;
cone_map.plot_skip      = 100   ;
cone_map.display_every  = 500   ;
cone_map.ID             = 0     ;
cone_map.max_time       = 200000;
cone_map.N_fast         = 1     ;

cone_map.initX = initialize_X( cone_map.M0, cone_map.M1, ...
                               cone_map.N_colors, cone_map.SS, ...
                               cone_map.cone_params.replusion_radii, ...
                               1, 1) ;
if ~isfield(cone_map.initX,'connections')
    cone_map.init.connections = zeros(N_GC,1) ;
end


% sanity check: compare cone_map.make_STA_W with make_LL
mLL = max(cone_map.LL(:)) ;
mx = 63 ;
my = 32 ;
mc = 1  ;
tX = flip_LL( cone_map.initX , [mx my mc] , cone_map , [1 1] ) ;
fprintf('\nLL and flip_ll: %f,%f, at x%d,y%d,c%d\n',mLL,tX.ll,mx,my,mc) ;
range_x = mx+(-4:5) ;
range_y = my+(-4:5) ;
test = zeros(numel(range_x),numel(range_y)) ;
for iii=range_x
    for jjj=range_y
        tX = flip_LL( cone_map.initX , [iii jjj 1] , cone_map , [1 1] ) ;
        test(iii,jjj) = tX.ll ;
    end
end

fprintf('\n')
disp(test(range_x,range_y))
disp(cone_map.LL(range_x,range_y,1))
disp( test(range_x,range_y) - cone_map.LL(range_x,range_y,1) )
fprintf('\n')
                           
end

function [sparse_struct, LL] = make_sparse_struct(cone_map,STA,WW,gaus_boxed)

M0 = cone_map.M0 ;
M1 = cone_map.M1 ;
SS = cone_map.SS ;
support = cone_map.cone_params.support_radius ;
colors  = cone_map.cone_params.colors ;

LL = zeros(M0*SS,M1*SS,cone_map.N_colors) ;

supersamples = 1/(2*SS):1/SS:1 ;
gs = cell(SS) ;
sparse_struct = cell( M0*SS, M1*SS, cone_map.N_colors ) ;

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
    gcLL = zeros(M0*SS,M1*SS,cone_map.N_colors) ;
    for ii=1:SS
        for jj=1:SS    
            CC = zeros(M0*M1,cone_map.N_colors) ;
            for color=1:cone_map.N_colors
                CCC = conv2( squeeze(STA(color,:,:,gc)), gs{ii,jj} ) ;
                CCC = CCC(support+1:M0+support,support+1:M1+support) ;
                CC(:,color) = CCC(:) ;
            end
            C = 0.5 * cone_map.quad_factors(gc) * (CC * colors').^2 / WW(ii,jj) ;
            C = max(0,C+0.5*cone_map.N_cones_terms(gc)) ;
            gcLL( ii:SS:M0*SS, jj:SS:M1*SS, :) = ...
                gcLL( ii:SS:M0*SS, jj:SS:M1*SS, :) + reshape(C,[M0 M1 3]) ;
        end
    end
    [x,yc] = find( gcLL>0 ) ;
    y = 1+mod(yc-1,M1*SS) ;
    c = ceil( yc/(M1*SS) ) ;
    for i=1:numel(x)
        sparse_struct{x(i),y(i),c(i)} = int16([sparse_struct{x(i),y(i),c(i)} gc]) ;
    end
    LL = LL + gcLL ;
    fprintf(' %d',gc)
end
end


function filter = make_filter_new(M0,M1,i,j,gaus_boxed, support)
filter = zeros(M0,M1) ;
[g,t,r,b,l] = filter_bounds( i, j, M0, M1, gaus_boxed, support) ;
filter(t:b,l:r) = g ;   
% filter is inverted; doesn't matter for the dot product calculation though
% filter = filter(end:-1:1,end:-1:1) ;  % uncomment to uninvert
end