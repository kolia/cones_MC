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
gaus_in_box = gaus_in_a_box( cone_params.sigma ) ;

% prior_cov   = cone_params.stimulus_variance^2*N_GC/sum(STA_norm.^2) ;
% prior_cov   = cone_params.stimulus_variance^2*(N_GC-1)/sum(STA_norm.^2) ;
prior_cov   = cone_params.stimulus_variance^2*(N_GC-1)/sum(STA_norm.^2) ;

cone_map.N_cones_term = sum( log( prior_cov) - log(cell_consts(:)+prior_cov) ) ;
cone_map.quad_factor  = N_spikes.^2 ./ (cell_consts+prior_cov) ;
% why was quad_factor changed?

% cone_map.A_factor  = N_spikes .* cone_params.stimulus_variance ./ ...
%                      (cell_consts+prior_cov) ;

cone_map.cov_factor = cell_consts+prior_cov ;

% % stereotyped cone receptive field
% s = cone_params.sigma ;
% R = cone_params.support_radius * SS ;
% cone_RF = exp(-0.5 * ((-R:R)/(SS*s)).^2)' * exp(-0.5 * ((-R:R)/(SS*s)).^2) ;
% cone_RF = cone_RF / norm(cone_RF(:)) ; clear s



% sparse int matrix, with number of out-of-border adjacencies
cone_map.outofbounds = sparse([],[],[],M0*SS,M1*SS,2*(M0+M1)*SS) ;
cone_map.outofbounds(:,[1 M1*SS]) = 1 ;
cone_map.outofbounds([1 M0*SS],:) = cone_map.outofbounds([1 M0*SS],:) + 1 ;




%% SETUP for Log-LIKELIHOOD calculations

% try load('STA_W')  % calculating STA_W takes a while, check for STA_W.mat
% catch e
% 
%     full_filters = cell(SS,SS) ;
% %     for i=
%     
%     %zeros(NROI*N_colors,N_GC) ;
%     factor = 0.5 ;
%     STA_W = sparse([],[],[],NROI*N_colors,N_GC,NROI*N_colors*N_GC*factor+10);
%     fprintf('Calculating STA_W...\n')
%     t = cputime ;
%     for ii=1:NROI
%         if mod(ii,100) == 0
%             fprintf('%d of %d in %.0f,  ',ii,NROI,cputime-t)
%         end
%         i   = cone_map.ROI(ii,1) ;
%         j   = cone_map.ROI(ii,2) ;
%         [index,filt]  = make_filter(M0,M1,i,j,cone_params.support_radius,gaus_in_box) ;
% 
%         for c=1:N_colors
%             filter  = kron(cone_params.colors(c,:),filt') ;
%             full_sta_w_ii = filter * STA([index index+M0*M1 index+2*M0*M1],:) ;
%             inds = abs(full_sta_w_ii) > quantile(abs(full_sta_w_ii), factor) ;
%             STA_W((c-1)*NROI + ii , inds) = full_sta_w_ii(inds) ;
%         end
%     end   % repmat( (N_spikes ./ cell_consts)' ,NROI*N_colors,1) ;
%     STA_W = STA_W ./ cone_params.stimulus_variance ; 
% 
% save('STA_W','STA_W')
% end
% cone_map.NICE           = NICE ;

cone_map.M0             = M0 ;
cone_map.M1             = M1 ;
cone_map.N_colors       = N_colors ;
cone_map.N_GC           = N_GC ;
cone_map.cell_consts    = cell_consts ;
cone_map.prior_cov      = prior_cov ;
cone_map.NROI           = NROI ;
cone_map.N_spikes       = N_spikes ;
cone_map.SS             = SS ;

coneConv = zeros( 2*R+SS , 2*R+SS , SS , SS ) ;
WW = zeros(SS,SS) ;

f = 1/(2*SS):1/SS:2*R/SS+1 ;
for xx=1:2*R+SS
    x = f(xx) ;
    for yy=1:2*R+SS
        y = f(yy) ;
        
        a = make_filter(4*R/SS+1,4*R/SS+1,x+R/SS,y+R/SS,cone_params.support_radius,gaus_in_box) ;
        
        for ss=1:SS
            s = (ss-0.5)/SS ;
            for tt=1:SS
                t = (tt-0.5)/SS ;
                
                b = make_filter(4*R/SS+1,4*R/SS+1,2*R/SS+s,2*R/SS+t,...
                    cone_params.support_radius,gaus_in_box) ;
                
                coneConv(xx,yy,ss,tt) = dot(a(:),b(:)) ;
            end
        end
        if (xx<=SS) && (yy<=SS)
            WW(xx,yy) = dot(a(:),a(:)) ;
        end
    end
end

cone_map.make_STA_W = memoized_STA_W(M0,M1,SS,cone_map.ROI,...
                                     cone_params.support_radius,cone_params.fudge,gaus_in_box) ;
                                 
cone_map.LL = make_LL(cone_map,STA,WW,gaus_in_box) ;

IC = inv(cone_params.colors) ;
QC = reshape( reshape(cone_map.LL,[],3) * IC' , size(cone_map.LL) ) ;
QC = QC - min(QC(:)) ;
cone_map.NICE = QC ./ max(QC(:)) ;

imagesc( cone_map.NICE )

% % test cone_map.make_STA_W against make_LL
% mLL = max(cone_map.LL(:)) ;
% [mk,mc] = find( reshape( cone_map.LL, NROI, N_colors) == mLL ) ;
% mx = mod( mk-1, M0*SS ) + 1 ;
% my = ceil( mk/(M0*SS) ) ;
% sta_w = cone_map.make_STA_W( mk, mc, reshape(STA,[],N_GC), cone_params.colors ) ;
% fprintf('\nLL and sta_w ll: %f,%f, %f, %f\n',mLL,cone_map.LL(mk+(mc-1)*NROI),...
%     0.5 * sum( cone_map.quad_factor' .* (sta_w / mean(WW(:))) .* sta_w  ),...
%     0.5 * sum( cone_map.quad_factor' * ((sta_w / mean(WW(:))) .* sta_w )')) ;
% 
% test = zeros(10,10) ;
% % range_x = 185:232 ;
% % range_y = 457:504 ;
% % mx = 194 ;
% % my = 465 ;
% % range_x = mx-5:mx+5 ; %1:10 ;
% % range_y = my-5:my+5 ; %1:10 ;
% range_x = 1:10 ;
% range_y = 1:10 ;
% for iii=range_x
%     for jjj=range_y
%         sta_w = cone_map.make_STA_W( iii+M0*SS*(jjj-1), 1, ...
%                                      reshape(STA,[],N_GC), cone_params.colors ) ;
%         test(iii,jjj) = 0.5 * sum( cone_map.quad_factor' .* ...
%                                   (sta_w / mean(WW(:))) .* sta_w ) ;
%     end
% end
% 
% disp(test(range_x,range_y))
% disp(cone_map.LL(range_x,range_y,1))

cone_map.R              = R ;
cone_map.coneConv       = coneConv ;
% cone_map.sumLconst      = sum(log(2*pi*cell_consts)) ;
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
                               cone_map.cone_params.replusion_radii,...
                               1, 1) ;
                           
end


% function LL = make_LL(M0,M1,N_colors,SS,STA,support,gaus_boxed)
function LL = make_LL(cone_map,STA,WW,gaus_boxed)

M0 = cone_map.M0 ;
M1 = cone_map.M1 ;
SS = cone_map.SS ;
N_colors = cone_map.N_colors ;
support = cone_map.cone_params.support_radius ;
colors  = cone_map.cone_params.colors ;

% supersamples = 1-1/(2*SS):-1/SS:0 ;
supersamples = 1/(2*SS):1/SS:1 ;
LL = zeros(M0*SS,M1*SS,N_colors) ;
for ii=1:SS
    for jj=1:SS
        i = supersamples(ii) ;
        j = supersamples(jj) ;
        ox  = 1+(floor(i-support):floor(i+support)) ;
        oy  = 1+(floor(j-support):floor(j+support)) ;
        x   = repmat(ox(:),numel(oy),1) ;
        y   = reshape( repmat(oy,numel(ox),1) , [] , 1 ) ;
        g   = reshape( gaus_boxed(i-x,j-y), [numel(ox) numel(oy)]) ;
        g   = g(end:-1:1,end:-1:1) ;
        for gc=1:cone_map.N_GC
            CC = zeros(M0*M1,N_colors) ;
            for color=1:N_colors
                CCC = conv2( STA(:,:,color,gc), g ) ;
                CCC = CCC(support+1:M0+support,support+1:M1+support) ;
                CC(:,color) = CCC(:) ;
            end
%             C = 0.5 * cone_map.quad_factor(gc) * (C * colors').^2 / WW(SS-ii+1,SS-jj+1) ;
            C = 0.5 * cone_map.quad_factor(gc) * (CC * colors').^2 / WW(ii,jj) ;
            LL( ii:SS:M0*SS, jj:SS:M1*SS, :) = ...
                LL( ii:SS:M0*SS, jj:SS:M1*SS, :) + reshape(C,[M0 M1 3]) ;
%             LL( SS-ii+1:SS:M0*SS, SS-jj+1:SS:M1*SS, :) = ...
%                 LL( SS-ii+1:SS:M0*SS, SS-jj+1:SS:M1*SS, :) + reshape(C,[M0 M1 3]) ;
        end
        fprintf('(%f, %f) ',ii,jj)
    end
end
end


function make_sta_w = memoized_STA_W(M0,M1,SS,ROI,support,fudge,gaus_boxed)
    
    NROI = M0*M1*SS*SS ;
    sta_w_memo = cell(NROI,3) ;
    usage = sparse([],[],[],NROI,3) ;
    filters = cell(SS,SS) ;
    
    make_sta_w = @make_STA_W ;
    function sta_w = make_STA_W( k, c, STA, colors )
        if isempty(sta_w_memo{k,c})
            i   = ROI(k,1) ;
            j   = ROI(k,2) ;

            ox    = 1+(floor(i-support):floor(i+support)) ;
            oy    = 1+(floor(j-support):floor(j+support)) ;
            ix    = (ox>=1) .* (ox<=M0) ;
            iy    = (oy>=1) .* (oy<=M1) ;
%             ox    = ox(ix>0) ;
%             oy    = oy(iy>0) ;
            ix    = repmat(ix(:),numel(iy),1) ;
            iy    = reshape( repmat(iy,numel(iy),1) , [] , 1 ) ;
            x     = repmat(ox(:),numel(oy),1) ;
            y     = reshape( repmat(oy,numel(ox),1) , [] , 1 ) ;

            ii = (i-floor(i))*SS+0.5 ;
            jj = (j-floor(j))*SS+0.5 ;
%             within = (i-support>=1) && (i+support<=M0) && (j-support>=1) && (j+support<=M1) ;
%             if within
%                 [ii jj]
%                 filters{ii,jj}
%             else
%                 [i j support]
%             end
            if ~isempty( filters{ii,jj} )
                filter = filters{ii,jj}.filter ;
            else
                filter = gaus_boxed( i-x, j-y ) ;
                filters{ii,jj}.filter = filter ;
            end

            index = (x+(y-1)*M0)' ;
            index = index((ix.*iy)>0) ;
            filter  = kron(colors(c,:),filter((ix.*iy)>0)') ;
            sta_w = (filter * STA([index index+M0*M1 index+2*M0*M1],:)) * fudge ;
            
            usage(k,c) = usage(k,c) + 1 ;
        else
            sta_w = sta_w_memo{k,c} ;
        end
        
        % purge least used half of sta_w_memo if more than 100000 values
        if sum(usage>0)>100000
            keep  = full( median(usage(:)) ) ;
            purge = find(usage<keep) ;
            sta_w_memo(purge) = cell(numel(purge),1) ;
            usage(purge) = 0 ;
        end
    end
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