function R = overlap_relation( X , Y )
% Calculate the binary relation R, a boolean square matrix, whose entry i,j
% is true iff cone i in configuration X is within the exclusion disk of
% cone j in configuration Y.

R       = sparse([],[],[],numel(X.state),numel(X.state),100) ;

[x,y]   = find(X.state) ;
xx      = find(X.state) ;

for i=1:length(x)
    [dum,indices] = place_mask( Y.M0 , Y.M1 , x(i) , y(i) , Y.masks{1,1}.exclusion ) ;
    R(xx(i),indices( Y.state(indices)>0 )) = true ;
end

end