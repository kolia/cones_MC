function X = make_contacts( X , x , y , id )

for d=1:4
    shift    = X.masks.cardinal{d} ;
    [dummy,inds] = place_mask( X.M0 , X.M1 , x , y , shift ) ;
    contacts    = X.id(inds) ;
    contacts    = contacts(contacts>0) ;
    
    X.contact{d}( id , contacts )             = true ;
    X.contact{1+mod(d+1,4)}( contacts , id )  = true ;
    X.contact{d}( id , id )                   = true ;
end

end