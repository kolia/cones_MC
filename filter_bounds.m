function [filter,t,r,b,l] = filter_bounds( i, j, M0, M1, gaus_boxed, support)

filter = gaus_boxed( i, j ) ;

b = min(M0,floor(i+support)+1) ;
r = min(M1,floor(j+support)+1) ;
t = max( 1,floor(i-support)+1) ;
l = max( 1,floor(j-support)+1) ;

if t>floor(i-support)+1
    filter = filter(t-floor(i-support):end,:) ;
end

if b<floor(i+support)+1
    filter = filter(1:end-(floor(i+support)+1-b),:) ;
end

if l>floor(j-support)+1
    filter = filter(:,l-floor(j-support):end) ;
end

if r<floor(j+support)+1
    filter = filter(:,1:end-(floor(j+support)+1-r)) ;
end

% [n,m] = size(filter) ;
% if n ~= b-t+1  || m ~= r-l+1
%     'asdfasdf'
% end

end