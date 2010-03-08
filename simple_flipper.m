function X = simple_flipper( X , flip_inds , ll )
% Simply flip some bits in X.state and recalculate log-likelihood from
% scratch.

flips = zeros(1,length(X.state)) ;
flips(flip_inds) = 1 ;
X.state = xor( X.state , flips ) ;
X.ll = ll(X.state) ;

end