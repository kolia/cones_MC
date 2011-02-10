function [sta,samples] = posterior_cone_map( cone_params , state , dX )

[N,M]   = size(state) ;
r = cell(4,1) ;
for c=1:3
    r{c}  = zeros(N,M) ;
end
r{4} = 0 ;

r = fold_states(state,dX,r,@accumulate_counts) ;

samples = r{4} ;

% gauscdf = gaus_in_a_box( cone_params.sigma ) ;

SS  = cone_params.supersample ;
R   = ceil( cone_params.support_radius * SS ) ;
x  = -R:1/SS:R ;
nn = length(x) ;
n  = floor(length(x)/2) ;
x  = repmat(x,[length(x) 1]) ;
xx = x' ;
gausfilter = mvnpdf([xx(:) x(:)],[],cone_params.sigma.*[1 1]) ;
gausfilter = reshape( gausfilter , [nn nn] ) ;
% gausfilter = gausfilter/sum(gausfilter(:)) ;

sta = zeros(N,M,3) ;
for c=1:3
    contrib_c = conv2(gausfilter,r{c}/r{4}) ;
    contrib_c = contrib_c(n+1:end-n,n+1:end-n) ;
    for cc=1:3
        sta(:,:,cc) = sta(:,:,cc) + contrib_c .* cone_params.colors(c,cc) ;
    end
end

end


function r = accumulate_counts( r , state )

% for each color, add counts for cones of that color.
for c=1:3
    r{c}(state == c) = r{c}(state == c) + 1 ;
end

% update number of samples seen so far
r{4} = r{4} + 1 ;

end