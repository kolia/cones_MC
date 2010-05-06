function [X,swaps] = swap_closure( X , otherX , XLL , OLL )
% Combine two systems into a single system, summing their log-likelihoods.
% This is in preparation for swapping parts of the two systems'
% configurations. For this, groups of cells that must be swapped together
% are calculated.

% calculate overlaps of X cones on otherX exclusion disks
R = overlap_relation( X , otherX ) ;

% symmetrize relation and calculate transitive closure
R = [ zeros(X.maxcones) R ; zeros(X.maxcones) R' ] ;
R = transitive_closure(R,1:X.maxcones,1) ;

% get equivalence classes / connected components
[classes,sizes] = equivalence_classes(R) ;

classes = classes( sizes>2 ) ;

swaps = cell( length(classes), 1 ) ;

for i=1:length(swaps)
    YX = X ;
        YO = otherX ;
        swaps{ns}.X     = YX ;
        swaps{ns}.with  = YO ;
        swaps{ns}.ll    = YX.ll + YO.ll ;
        swaps{ns}.state = YX.state ;
        swaps{ns}.forward_prob = 1 ;  % all swaps have the same prob anyway
end
        
X.X     = X ;
X.with  = otherX ;
X.ll    = X.X.ll + X.with.ll ;
X.state = X.X.state ;

if length(swaps) < 2
    swaps{1}.type = 'whole' ;
end

end