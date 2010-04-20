function [mask,indices] = place_mask( M0 , M1 , x , y , mask )

mask    = [ mask(:,1)+x   mask(:,2)+y ] ;
lines   = mask(:,1)>0  &  mask(:,1)<=M0  &  mask(:,2)>0  &  mask(:,2)<=M1 ;
mask    = mask(lines,:) ;

if nargout>1
    indices = mask(:,1) + (mask(:,2)-1)*M0 ;
end

end