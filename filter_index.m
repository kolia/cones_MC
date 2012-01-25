function [filter,index] = filter_index( i, j, M0, M1, gaus_boxed, support)

filter = gaus_boxed( i, j ) ;

ox    = 1+(floor(i-support):floor(i+support)) ;
oy    = 1+(floor(j-support):floor(j+support)) ;

ix     = (ox>=1) .* (ox<=M0) ;
iy     = (oy>=1) .* (oy<=M1) ;
ix     = repmat(ix(:),numel(iy),1) ;
iy     = reshape( repmat(iy,numel(iy),1) , [] , 1 ) ;

x      = repmat(ox(:),numel(oy),1) ;
y      = reshape( repmat(oy,numel(ox),1) , [] , 1 ) ;

index  = (x+(y-1)*M0)' ;
index  = index((ix.*iy)>0) ;
filter = filter((ix.*iy)>0)' ;

end