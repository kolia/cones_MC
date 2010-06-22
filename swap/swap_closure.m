function XS = swap_closure( X , otherX , PROB )
% Combine two systems into a single system, summing their log-likelihoods.
% This is in preparation for swapping parts of the two systems'
% configurations. For this, groups of cells that must be swapped together
% are calculated.

% calculate overlaps of X cones on otherX exclusion disks
R = overlap_relation( otherX , X ) ;

% symmetrize relation and calculate transitive closure
N = numel(X.state) ;

i = find( otherX.state) ;
EO= logical(sparse(i,i,ones(length(i),1),N,N,3*N)) ;

j = find( X.state) ;
EX= logical(sparse(j,j,ones(length(j),1),N,N,3*N)) ;

R = logical( [ EO R ; R' EX ] ) ;
R = transitive_closure(R,[i ; N+j],0) ;

clear i j

% get equivalence classes / connected components
[classes,sizes] = equivalence_classes(R) ;
classes = classes( sizes>2 ) ;
classes = [{[]} ; classes] ;
XS      = cell( length(classes), 1 ) ;

% fprintf('\n#swaps: %d',length(XS))

keep = [] ;
for i=1:length(XS)
    XX = X ;
    OX = otherX ;
   
    Oclass   = classes{i}(classes{i}<=N) ;
    Xclass   = classes{i}(classes{i} >N) - N ;
    
%     fprintf('\t %d,%d',numel(Xclass),numel(Oclass)) ;
        
    xcx      = 1+mod(Xclass-1,X.M0) ;
    xcy      = 1+floor((Xclass-1)/X.M0) ;
    xcc      = XX.state(Xclass) ;
    
    ocx      = 1+mod(Oclass-1,X.M0) ;
    ocy      = 1+floor((Oclass-1)/X.M0) ;
    occ      = OX.state(Oclass) ;
    
    if X.maxcones >= X.N_cones      + numel(Oclass)    && ...
       X.maxcones >= otherX.N_cones + numel(Xclass)

            for k=1:length(Xclass)
                XX = flip_LL( XX , [xcx(k) xcy(k) 0] , PROB ) ;
            end
            for k=1:length(Oclass)
                OX = flip_LL( OX , [ocx(k) ocy(k) 0] , PROB ) ;
            end
            for k=1:length(Oclass)
                XX = flip_LL( XX , [ocx(k) ocy(k) occ(k)] , PROB ) ;
            end
            for k=1:length(Xclass)
                OX = flip_LL( OX , [xcx(k) xcy(k) xcc(k)] , PROB ) ;
            end
            
            XS{i}.X     = XX ;
            XS{i}.with  = OX ;
            XS{i}.ll    = XX.ll + OX.ll ;
            
        keep = [keep i] ;
    end
end
XS = XS(keep) ;

N = length(XS) ;
for i=1:N
    XS{i}.forward_prob = 1/N ;
end

end