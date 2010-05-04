function check_X( LL , Y )

% Y = X ;

% [x,y] = find( X.id ) ;

% % check X.contact
% for d=1:4
%     X.contact{d} = X.contact{d} * 0 ;
%     
%     for i=1:length(x)
%         X = make_contacts( X , x(i) , y(i) , X.id(x(i),y(i)) , LL ) ;
%     end    
%     
%     dX = find( xor( Y.contact{d} , X.contact{d}) ) ;
%     if dX
%         Z.contact{d}
%         X.contact{d}
%         Y.contact{d}
%         dX
%         d        
%     end
% end


% for i=1:length(x)
%     X = make_contacts( X , x(i) , y(i) , X.id(x(i),y(i)) , LL ) ;
% % check exclusion
% %     [~,indices] = place_mask( X.M0 , X.M1 , x(i) , y(i) , X.masks.exclusion ) ;
% %     
% %     overlap = find( X.state(indices)>0 , 1) ;
% %     if ~isempty( overlap )  &&   indices(overlap) ~= x(i) + X.M0*(y(i)-1)
% %         [x(i) y(i)]
% %     end
% end


% % check X.shift_dLL
% for d=1:4
%     f =  find( X.shift_dLL{d} ~= Y.shift_dLL{d} ) ;
%     if ~isempty(f)
%         'blah'
%     end
% end


[x,y] = find( Y.id ) ;

% calculate Y.ll directly
inds = x + (y-1)*Y.M0 ;
inds = inds + Y.M0*Y.M1*( Y.state(inds)-1 ) ;
Yll  = sum( LL(inds) ) - Y.N_cones_factor * length(x) ;

if abs( Y.ll - Yll ) > 1
    'aha!'
end

end