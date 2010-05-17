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

% check that localLL is consistent with LL(inds)
[Xx,Xy,Xc] = find(Y.state) ;
[x,y,v]    = find(Y.id   ) ;
inds   = Xx + Y.M0*(Xy-1) + Y.M0*Y.M1*(Xc-1) ;
if ~isempty(find( Y.localLL(v) - LL(inds) , 1))
   'ha!' 
end

% check if N_cones if consistent with id, state and taken_ids
counts = [numel(find(Y.localLL)) numel(find(Y.id)) Y.N_cones numel(find(Y.state)) numel(find(Y.taken_ids))] ;
if ~isempty(find(diff(counts)))
    counts
    'N_cones is off'
end

% calculate Y.ll directly
Yll  = sum( Y.localLL ) - Y.N_cones_factor * Y.N_cones ;
if abs( Y.ll - Yll ) > 1
    'aha!'
end

end