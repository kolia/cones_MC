function R = overlap_relation( X , Y )

R    = false(X.maxcones) ;
mask = Y.masks.exclusion ;

[x,y,v]= find(X.id) ;

for i=1:length(x)
    [~,indices] = place_mask( X.M0 , X.M1 , x(i) , y(i) , X.masks.exclusion ) ;
    
    ids = Y.id( Y.id(indices)>0 ) ;
    
    R(v(i),ids) = true ;
end

end