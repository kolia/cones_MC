function X = make_contacts( X , x , y , id , LL )

if nargin < 5
    id = X.id(x,y) ;
end

for d=1:4
    shift    = X.masks.cardinal{d} ;
    [~,inds] = place_mask( X.M0 , X.M1 , x , y , shift ) ;
    contacts    = X.id(inds) ;
    contacts    = contacts(contacts>0) ;
    
    X.contact{d}( id , contacts )             = true ;
    X.contact{1+mod(d+1,4)}( contacts , id )  = true ;
    
    % calculate local dLL
    X.shift_dLL{d}(id) = local_dLL(LL,x,y,shift,X.M0,X.M1) ;
end

end