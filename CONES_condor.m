function ticket = CONES_condor( GC_stas , cone_params , type , min_overlap , s )
% ticket = CONES_condor( GC_stas , cone_params , type , min_overlap , s )
% launches condor cone finding job

if nargin<5 ,  s = 12 ; end
if nargin<4 ,  min_overlap  = 6 ; end
if nargin<3 ,  type         = 'MCMC' ; end

% make sure agricola is present on the current path
[dummy,o] = unix('ls -l | grep ''agricola$''') ;
if isempty(o)
    agricola_folder = which('sow') ;
    if isempty(agricola_folder)
        error('agricola folder not currently in matlab path!') ;
    end
    agricola_folder = agricola_folder(1:end-6) ;
    unix(sprintf('cp -r %s .',agricola_folder)) ;
end

% cut problem up into regions
[M0,M1,N_color] = size(GC_stas(1).spatial) ;
cutX = cut_interval( M0 , s , min_overlap ) ;
cutY = cut_interval( M1 , s , min_overlap ) ;

nX   = numel(cutX) ;
nY   = numel(cutY) ;

cone_maps = cell(nX,nY) ;
for i=1:nX
    for j=1:nY
        cone_maps{i,j}{1}.cut.X           = cutX ;
        cone_maps{i,j}{1}.cut.Y           = cutY ;
        cone_maps{i,j}{1}.cut.position    = { cutX(i)+(0:s-1) ; cutY(j)+(0:s-1) ; 1:3 } ;
        cone_maps{i,j}{1}.cut.total_sizes = [M0 M1 N_color] ;
        cone_maps{i,j}{1}.ROI             = zeros(M0,M1,N_color) ;
        cone_maps{i,j}{1}.ROI(cutX(i)+(0:s-1),cutY(j)+(0:s-1),1:3) = 1 ;
    end
end
cone_maps = cone_maps(:) ;

% launch cluster of jobs
switch type
    case 'MCMC'
        ticket = sow(sprintf('cones_%d_%d_MC_%d',s,min_overlap,cone_params.supersample) , @(conemap)MCMC_cones(GC_stas,cone_params,conemap) , cone_maps) ;
        return
    case 'greedy'
        ticket = sow(sprintf('cones_%d_%d_GREEDY_%d',s,min_overlap,cone_params.supersample) , @(conemap)greedy_cones(GC_stas,cone_params,conemap) , cone_maps) ;
        return
    otherwise
        error('CONES_condor( GC_stas , cone_params , type ) : type must be ''MCMC'' or ''greedy''')
end

end

function cuts = cut_interval( N , s , min_overlap )
n       = ceil(  (N-1 - min_overlap) / (s-min_overlap) ) ;
cuts    = floor(1:(N-s)/(n-1):N-s+1) ;
end