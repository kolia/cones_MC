function X = make_contacts( X , xyc )

x = xyc(1) ;
y = xyc(2) ;
c = xyc(3) ;

xy           = x+(y-1)*X.M0 ;

for d=1:2
    shift        = X.masks{1,1}.cardinal{d} ;
    [dummy,inds] = place_mask( X.M0 , X.M1 , x , y , shift ) ;
    contacts     = inds(X.state(inds)>0) ;    
    
    X.contact{d}( xy , contacts )             = true ;

    shift        = X.masks{1,1}.cardinal{d+2} ;
    [dummy,inds] = place_mask( X.M0 , X.M1 , x , y , shift ) ;
    contacts     = inds(X.state(inds)>0) ;
    X.contact{d}( contacts , xy )             = true ;
    
%     X.contact{1+mod(d+1,4)}( contacts , xy )  = true ;
%     X.contact{d}( xy , xy )                   = true ;

%     if ~X.state(contacts)
%         X.state
%         contacts
%     end

end

% for d=1:2
%     if nnz(X.contact{d} - X.contact{1+mod(d+1,4)}')>0
%         'oups1'
%     end
% end

end