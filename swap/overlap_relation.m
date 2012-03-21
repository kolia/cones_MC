function R = overlap_relation( X , Y )
% Calculate the binary relation R, a boolean square matrix, whose entry i,j
% is true iff cone i in configuration X is within the exclusion disk of
% cone j in configuration Y.

[x,y]   = find(X.state) ;
xx      = find(X.state) ;

x_inds = zeros(numel(x),1) ;
y_inds = zeros(numel(x),1) ;
ind    = 0 ;

for i=1:length(x)
    [dum,indices] = place_mask( Y.M0 , Y.M1 , x(i) , y(i) , Y.masks{1,1}.exclusion ) ;
    overlaps = find( Y.state(indices)>0 ) ;
    for j=1:numel(overlaps)
        ind = ind+1 ;
        x_inds(ind) = xx(i) ;
        y_inds(ind) = indices(overlaps(j)) ;
    end
end

x_inds = x_inds(1:ind) ;
y_inds = y_inds(1:ind) ;

R = logical( sparse(x_inds,y_inds,ones(ind,1),numel(X.state),numel(X.state))) ;

end