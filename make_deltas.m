function deltas = make_deltas(m,M,r,N)

N = max(N-1,4) ;

deltas = log( 1:exp(r)/(N-1):1+exp(r) ) ;
deltas = m + deltas*(M-m)/deltas(end) ;
deltas = [deltas(1:2) (deltas(2)+deltas(3))/2 deltas(3:end-3) (deltas(end-1)+deltas(end-2))/2 deltas(end-1:end)] ;

deltas = deltas(end:-1:1) ;

end