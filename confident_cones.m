function cc = confident_cones( X , dX , PROB , selector )
% [sta,samples] = denoised_sta( X , dX , PROB , truecolor )
% Integrate over MCMC sample path starting at X with moves dX, to produce
% denoised STAs for all ganglion cells.  PROB contains parameters, as
% output by exact_LL_setup.m
% If the optional truecolor is true, then true RGB STAs are output,
% otherwise 


if nargin<4
    selector = @(n) (n>10000) && (mod(n,20) == 0) ;
end
    
N = PROB.M0 * PROB.SS ;
M = PROB.M1 * PROB.SS ;
r = cell(3,1) ;
r{1} = zeros(N,M,3) ;
r{2} = 0 ;
r{3} = 0 ;

r = fold_X(X,dX,PROB,[1 1],r,@(r,x)accumulate_stas(r,x,selector)) ;

cc = r{1}/r{2} ;

end


function r = accumulate_stas( r , X , selector )

% only accumulate when selector is true
if selector(r{3})
    for c=1:3
        r{1}(:,:,c) = r{1}(:,:,c) + (X.state==c) ;
    end
        
    % update number of samples accumulated so far
    r{2} = r{2} + 1 ;
end
% update number of samples replayed so far
r{3} = r{3} + 1 ;

if ~mod(r{3},1000)
    fprintf('.')
end

end