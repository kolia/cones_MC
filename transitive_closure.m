function R = transitive_closure( R , inds , sym)
% Calculate transitive closure of a binary relation
% represented by boolean matrix R.
% inds is a list of indices which are known to be necessary for this
% calculation. 2 cases:
% - transitive_closure is being called for the first time, inds should be
% the indices of rows and columns where R has at least one true entry.
% - transitive_closure has already been called, and the only modification
% to R has been to add some entries, then inds need only be the indices of
% those entries.

% non-symmetrical relation
if nargin<3 || ~sym
    for k=inds(:)'
        R( R(:,k) , R(k,:) ) = true ;
    end
    
% symmetrical relation
else    
    for k=inds(:)'
        i = R(:,k) ;
%         if ~isempty(find(i,1))
            R( i , i ) = true ;
%         end
    end
    
end


end