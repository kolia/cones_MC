function [sta,invww] = denoised_sta( X , dX , PROB , selector, signif , truecolor )
% [sta,samples] = denoised_sta( X , dX , PROB , truecolor )
% Integrate over MCMC sample path starting at X with moves dX, to produce
% denoised STAs for all ganglion cells.  PROB contains parameters, as
% output by exact_LL_setup.m
% If the optional truecolor is true, then true RGB STAs are output,
% otherwise 


if nargin<6
    truecolor = false ;
end

if nargin<5
    signif = 0 ;
end

if nargin<4
    selector = @(n) (n>10000) && (mod(n,20) == 0) ;
end
    
N = PROB.M0 * PROB.SS ;
M = PROB.M1 * PROB.SS ;
r = cell(5,1) ;
for c=1:3
    r{c} = cell(PROB.N_GC,2) ;
    for i=1:PROB.N_GC
        r{c}{1,i} = zeros(N,M) ;
        r{c}{2,i} = zeros(N,M) ;
    end
end
r{4} = 0 ;
r{5} = 0 ;
r{6} = 0 ;

% gauscdf = gaus_in_a_box( cone_params.sigma ) ;

if ~isfield(X,'invWW') && X.N_cones==0
    X.invWW = [] ;
    if isfield(X,'WW')
        X = rmfield(X,'WW') ;
    end
end

R  = ceil( PROB.cone_params.support_radius * PROB.SS ) ;
x  = -R:1/PROB.SS:R ;
nn = length(x) ;
n  = floor(length(x)/2) ;
x  = repmat(x,[length(x) 1]) ;
xx = x' ;
gausfilter = mvnpdf([xx(:) x(:)],[], PROB.cone_params.sigma.*[1 1] ) ;
gausfilter = reshape( gausfilter , [nn nn] ) ;
% gausfilter = gausfilter/sum(gausfilter(:)) ;

r = fold_X(X,dX,PROB,[1 1],r,@(r,x)accumulate_stas(r,x,PROB,gausfilter,n,selector,signif)) ;

sta = cell(2,PROB.N_GC) ;
for i=1:PROB.N_GC
    sta{1,i} = zeros(N,M,3) ;
    sta{2,i} = zeros(N,M,3) ;
    for c=1:3
        contrib_c = r{c}{1,i}/r{4} ;
        contrib_sq = r{c}{2,i}/r{4} ;
        if truecolor
            for cc=1:3
                sta{1,i}(:,:,cc) = sta{1,i}(:,:,cc) + ...
                    contrib_c .* PROB.cone_params.colors(c,cc) ;
                sta{2,i}(:,:,cc) = sta{2,i}(:,:,cc) + ...
                    contrib_sq .* PROB.cone_params.colors(c,cc) ;
            end
        else
            sta{1,i}(:,:,c) = sta{1,i}(:,:,c) + contrib_c ;
            sta{2,i}(:,:,c) = sta{2,i}(:,:,c) + contrib_sq ;
        end
    end
end

invww = r{6}/(r{4}*PROB.cov_factor) ;

end


function r = accumulate_stas( r , X , PROB , gausfilter, n, selector, signif )

% only accumulate when selector is true
if selector(r{5})

    NGC = size(PROB.STA_W,2) ;

    invWW = X.invWW ;
    invWW(abs(invWW)<abs(invWW(1,1))*1e-17) = 0 ;

    invWW = sparse(invWW) ;

    [x,y,c] = find(X.state) ;

    inds = x+PROB.M0*PROB.SS*(y-1)+PROB.M0*PROB.M1*PROB.SS*PROB.SS*(c-1) ;
    STA_W_state = PROB.STA_W( inds , : ) ;

    factor = repmat( PROB.cell_consts' ./ PROB.cov_factor' , [numel(x) 1] ) ;
    A   = (invWW * STA_W_state) .* factor ;

    STD = sqrt(diag(invWW) * (1./PROB.cov_factor')) ;
    B   = A ./ STD ;

    A( abs(B) < signif ) = 0 ;

    % for each GC, add sample contribution.
    for i=1:NGC
        for cc=1:3
            I = c == cc ;
            if sum(I)>0
                pixels = zeros(PROB.M0*PROB.SS,PROB.M1*PROB.SS) ;
                pixels(x(I)+PROB.M0*PROB.SS*(y(I)-1)) = A(I,i) ;
                pixels = conv2(gausfilter,pixels)   ;
                pixels = pixels(n+1:end-n,n+1:end-n);
                r{cc}{1,i} = r{cc}{1,i} + pixels    ;
                r{cc}{2,i} = r{cc}{2,i} + pixels.^2 ;
    %             r{cc}{1,i}(x(I)+PROB.M0*PROB.SS*(y(I)-1)) = ...
    %                     r{cc}{1,i}(x(I)+PROB.M0*PROB.SS*(y(I)-1)) + pixels ;
    %             r{cc}{2,i}(x(I)+PROB.M0*PROB.SS*(y(I)-1)) = ...
    %                     r{cc}{2,i}(x(I)+PROB.M0*PROB.SS*(y(I)-1)) + a.^2 ;
            end
        end
    end

    % accumulate mean of diag( invWW )
    r{6} = r{6} + mean(diag(invWW)) ;
    
    % update number of samples accumulated so far
    r{4} = r{4} + 1 ;
end
% update number of samples replayed so far
r{5} = r{5} + 1 ;
end