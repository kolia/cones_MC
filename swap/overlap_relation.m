function R = overlap_relation( X , Y )
% Calculate the binary relation R, a boolean square matrix, whose entry i,j
% is true iff cone i in configuration X is within the exclusion disk of
% cone j in configuration Y.

R    = false(X.maxcones) ;

[x,y,v]= find(X.id) ;

for i=1:length(x)
    [dum,indices] = place_mask( Y.M0 , Y.M1 , x(i) , y(i) , Y.masks.exclusion ) ;
    ids = Y.id( indices( Y.id(indices)>0 ) ) ;

    R(v(i),ids) = true ;
end

end