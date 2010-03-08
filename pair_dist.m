function D = pair_dist(x)

N = size(x,2) ;
D = zeros(N) ;
for i=1:N
    for j=i:N
        d = norm( x(:,i) - x(:,j) ) ;
        D(i,j) = d ;
        D(j,i) = d ;
    end
end

end