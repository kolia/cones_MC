function deltas = make_deltas(m,M,r,N)
% make a progression of N temperatures in range [m M] with param r
% for example try:  make_deltas(0.1, 1, 1, 20 ) ;

N = max(N,4) ;

deltas = log( 1:exp(r)/(N-1):1+exp(r) ) ;
deltas = m + deltas*(M-m)/deltas(end) ;

deltas = deltas(end:-1:1) ;

end