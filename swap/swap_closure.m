function [result_X,swaps] = swap_closure( X , otherX , XLL , OLL )
% Combine two systems into a single system, summing their log-likelihoods.
% This is in preparation for swapping parts of the two systems'
% configurations. For this, groups of cells that must be swapped together
% are calculated.

% calculate overlaps of X cones on otherX exclusion disks
R = overlap_relation( X , otherX ) ;

% symmetrize relation and calculate transitive closure
R = logical( [ eye(X.maxcones) R ; R' eye(X.maxcones) ] ) ;
R = transitive_closure(R,1:2*X.maxcones,1) ;

% get equivalence classes / connected components
[classes,sizes] = equivalence_classes(R) ;

classes = classes( sizes>2 ) ;

swaps = cell( length(classes), 1 ) ;


[Xx,Xy,Xc] = find(X.state) ;
[dum,dum,Xi]   = find(X.id) ;
XX         = sortrows([Xx Xy Xc Xi],4) ;
Xinds      = zeros(X.maxcones,1) ;
Xinds(XX(:,4)) = XX(:,1) + X.M0*(XX(:,2)-1) + X.M0*X.M1*(XX(:,3)-1) ;

[Ox,Oy,Oc] = find(otherX.state) ;
[dum,dum,Oi]   = find(otherX.id) ;
OO         = sortrows([Ox Oy Oc Oi],4) ;
Oinds      = zeros(X.maxcones,1) ;
Oinds(OO(:,4)) = OO(:,1) + otherX.M0*(OO(:,2)-1) + otherX.M0*otherX.M1*(OO(:,3)-1) ;

if length(swaps) > 1
    for i=1:length(swaps)
        Xclass   = classes{i}(classes{i}<=X.maxcones) ;
        Oclass   = classes{i}(classes{i} >X.maxcones) - X.maxcones ;
        
        swaps{i} = make_swap( X , otherX , XLL , OLL , Xinds , Oinds , Xclass , Oclass , 'patch' ) ;
    end
else
    Xclass   = find(     X.taken_ids) ;
    Oclass   = find(otherX.taken_ids) ;
    swaps{1} = make_swap( X , otherX , XLL , OLL , Xinds , Oinds , Xclass , Oclass , 'whole' ) ;
end

result_X.X     = X ;
result_X.with  = otherX ;
result_X.ll    = X.ll + otherX.ll ;
result_X.state = result_X.X.state ;

end


function swap = make_swap( X , otherX , XLL , OLL , Xinds , Oinds , Xclass , Oclass , type )

dN          = length(Xclass) - length(Oclass) ;

LL_X        = X.ll      + dN *      X.N_cones_factor ...
            - sum(      X.localLL(Xclass)) + sum( XLL(Oinds(Oclass))) ;
LL_O        = otherX.ll - dN * otherX.N_cones_factor ...
            - sum( otherX.localLL(Oclass)) + sum( OLL(Xinds(Xclass))) ;

swap.ll = LL_X + LL_O ;
swap.forward_prob = 1 ;  % all swaps have the same trial prob
swap.move = {type Xclass Oclass LL_X LL_O} ;

end