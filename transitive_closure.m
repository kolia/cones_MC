function R = transitive_closure( R , inds )
% Calculate transitive closure of a binary relation
% represented by boolean matrix R.
% inds is a list of indices which are known to be necessary for this
% calculation. 2 cases:
% - transitive_closure is being called for the first time, inds should be
% the indices of rows and columns where R has at least one true entry.
% - transitive_closure has already been called, and the only modification
% to R has been to add some entries, then inds need only be the indices of
% those entries.

for k=inds(:)'
    R( R(:,k) , R(k,:) ) = true ;
end

end