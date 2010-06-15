function X = make_contacts( X , x , y )

for d=1:4
    shift        = X.masks.cardinal{d} ;
    [dummy,inds] = place_mask( X.M0 , X.M1 , x , y , shift ) ;
    contacts     = inds(X.state(inds)>0) ;
    
    xy           = x+(y-1)*X.M0 ;
    
    X.contact{d}( xy , contacts )             = true ;
    X.contact{1+mod(d+1,4)}( contacts , xy )  = true ;
    X.contact{d}( xy , xy )                   = true ;

%     if ~X.state(contacts)
%         X.state
%         contacts
%     end

end

end