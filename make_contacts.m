function X = make_contacts( X , xyc )

x = xyc(1) ;
y = xyc(2) ;
c = xyc(3) ;

for d=1:4
    shift        = X.masks{1,1}.cardinal{d} ;
    [dummy,inds] = place_mask( X.M0 , X.M1 , x , y , shift ) ;
    contacts     = inds(X.state(inds)>0) ;
    
    xy           = x+(y-1)*X.M0 ;
    
    if d<=2
        X.contact{d}( xy , contacts )             = true ;
    else
        X.contact{1+mod(d+1,4)}( contacts , xy )  = true ;
    end
%     X.contact{1+mod(d+1,4)}( contacts , xy )  = true ;
%     X.contact{d}( xy , xy )                   = true ;

%     if ~X.state(contacts)
%         X.state
%         contacts
%     end

end

end