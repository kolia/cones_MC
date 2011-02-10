function X = downsample(X,n)

X = ds(X,n) ;
X = X' ;
X = ds(X,n) ;
X = X' ;

end


function X = ds(X,n)

N = size(X,2) ;
if N/n ~= ceil(N/n)
    error('downsampling must be for an exact multiple.')
end

down = zeros(N,N/n) ;
for i=1:N/n
    down(n*(i-1)+1:n*i,i) = 1 ;
end

X = X*down/n ;

end