function [sta,r] = denoised_sta( X , dX , PROB , signif , truecolor )
% [sta,samples] = denoised_sta( X , dX , PROB , truecolor )
% Integrate over MCMC sample path starting at X with moves dX, to produce
% denoised STAs for all ganglion cells.  PROB contains parameters, as
% output by exact_LL_setup.m
% If the optional truecolor is true, then true RGB STAs are output,
% otherwise 


if nargin<5
    truecolor = true ;
end

if nargin<4
    signif = 0 ;
end
    
N = PROB.M0 * PROB.SS ;
M = PROB.M1 * PROB.SS ;
r = cell(4,1) ;
for c=1:3
    r{c} = cell(PROB.N_GC,1) ;
    for i=1:PROB.N_GC
        r{c}{i} = zeros(N,M) ;
    end
end
r{4} = 0 ;

r = fold_X(X,dX,PROB,r,@(r,x)accumulate_stas(r,x,PROB,signif)) ;

% gauscdf = gaus_in_a_box( cone_params.sigma ) ;

R  = ceil( PROB.cone_params.support_radius * PROB.SS ) ;
x  = -R:1/PROB.SS:R ;
nn = length(x) ;
n  = floor(length(x)/2) ;
x  = repmat(x,[length(x) 1]) ;
xx = x' ;
gausfilter = mvnpdf([xx(:) x(:)],[], PROB.cone_params.sigma.*[1 1] ) ;
gausfilter = reshape( gausfilter , [nn nn] ) ;
% gausfilter = gausfilter/sum(gausfilter(:)) ;

sta = cell(PROB.N_GC,1) ;
for i=1:PROB.N_GC
    sta{i} = zeros(N,M,3) ;
    for c=1:3
        contrib_c = conv2(gausfilter,r{c}{i}/r{4}) ;
        contrib_c = contrib_c(n+1:end-n,n+1:end-n) ;
        if truecolor
            for cc=1:3
                sta{i}(:,:,cc) = sta{i}(:,:,cc) + ...
                    contrib_c .* PROB.cone_params.colors(c,cc) ;
            end
        else
            sta{i}(:,:,c) = sta{i}(:,:,c) + contrib_c ;
        end
    end
end

end


function r = accumulate_stas( r , X , PROB , signif )

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
        r{cc}{i}(x(I)+PROB.M0*PROB.SS*(y(I)-1)) = ...
                r{cc}{i}(x(I)+PROB.M0*PROB.SS*(y(I)-1)) + A(I,i) ;
    end
end

% update number of samples seen so far
r{4} = r{4} + 1 ;

end