function check_X( X , LL , Z )

Y = X ;

for d=1:4
    X.contact{d} = X.contact{d} * 0 ;
    
    [x,y] = find( X.id ) ;
    for i=1:length(x)
        X = make_contacts( X , x(i) , y(i) , X.id(x(i),y(i)) , LL ) ;
    end    
    
    dX = find( xor( Y.contact{d} , X.contact{d}) ) ;
    if dX
        Z.contact{d}
        X.contact{d}
        Y.contact{d}
        dX
        d        
    end
    
end

for i=1:length(x)
    X = make_contacts( X , x(i) , y(i) , X.id(x(i),y(i)) , LL ) ;
    
    [~,indices] = place_mask( X.M0 , X.M1 , x(i) , y(i) , X.masks.exclusion ) ;
    
    overlap = find( X.state(indices)>0 , 1) ;
    if ~isempty( overlap )  &&   indices(overlap) ~= x(i) + X.M0*(y(i)-1)
        [x(i) y(i)]
    end

end


end