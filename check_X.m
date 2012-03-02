function check_X( Y )

% global cone_map

X = Y ;
[x,y,c] = find( X.state ) ;

if numel(find(X.state>0)) ~= size(X.WW,1)
   disp('X.state and X.WW have inconsistent sizes.')
end

if size(X.sparse_STA_W_state) ~= size(X.WW,1)
   disp('X.STA_W_state and X.WW have inconsistent sizes.')
end


% % check X.contact
% for d=1:4
%     X.contact{d} = X.contact{d} * 0 ;
%     
%     for i=1:length(x)
%         X = make_contacts( X , [x(i) y(i) c(i)]) ;
%     end    
%     
%     dX = find( xor( Y.contact{d} , X.contact{d}) ) ;
%     if dX
% %         X.diff.added
% %         X.diff.deleted
%         X.contact{d}
%         Y.contact{d}
%         dX
%         d
%     end
% end

% check that X.contact only concerns existing cones
for d=1:2
    cones = logical( max( max( X.contact{d} ,[],1) , max( X.contact{d} ,[],2)') ) ;
%     if ~min(X.state(cones))
%         X.state
%         cones
%     end
%     if ~cones(X.state(:)>0)
%         cones
%         X.state
%     end
end

% % check that ll is consistent   % NEEDS cone_map TO BE GLOBAL
% if isfield(X,'contact') ,  X = rmfield(X,'contact') ; end
% X.N_cones   = 0  ;
% X.invWW     = [] ;
% X.state     = 0 * X.state ;
% for k=1:length(x)
%     X = flip_LL( X , [x(k) y(k) c(k)] , cone_map ) ;
% end
% 
% if abs(X.ll - Y.ll)>1e-8
%     X.ll - Y.ll
%     X.diff
%     [x y c]
%     'll is off'
% end

for i=1:length(x)
%     X = make_contacts( X , x(i) , y(i) , X.id(x(i),y(i)) , LL ) ;
% check exclusion
    [dummy,indices] = place_mask( X.M0 , X.M1 , x(i) , y(i) , X.masks{1,1}.exclusion ) ;
    
    overlap = find( X.state(indices)>0 , 1) ;
    if ~isempty( overlap )  &&   indices(overlap) ~= x(i) + X.M0*(y(i)-1)
        [x(i) y(i) c(i)]
    end
end


% % check that localLL is consistent with LL(inds)
% [Xx,Xy,Xc] = find(Y.state) ;
% [x,y,v]    = find(Y.id   ) ;
% inds   = Xx + Y.M0*(Xy-1) + Y.M0*Y.M1*(Xc-1) ;
% if ~isempty(find( Y.localLL(v) - LL(inds) , 1))
%    'ha!' 
% end

% check if N_cones if consistent with id, state and taken_ids
% counts = [Y.N_cones numel(find(Y.id)) numel(find(Y.state)) ...
%           numel(find(Y.taken_ids)) size(Y.invWW,1)] ;
counts = [numel(find(Y.state)) size(Y.WW,1)] ;
if ~isempty(find(diff(counts)))
    counts
    'N_cones is off'
end


end