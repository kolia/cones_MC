function XS = swap_closure( X , T, otherX , oT, PROB )
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

% classes = [{[]} ; {[i' N+j']}] ;
classes = [{[]} ; classes] ;

XS      = cell( length(classes), 1 ) ;

% fprintf('\n#swaps: %d',length(XS))

keep = [] ;
for m=1:length(XS)
    if m>1
        XX       = struct ;
        XX.state = X.state ;
        XX.invWW = X.invWW ;
        XX.N_cones=X.N_cones ;
        XX.ll    = X.ll ;
        XX.diff  = X.diff ;
        XX.beta  = X.beta ;
        XX.delta = X.delta ;
        
        OX       = struct ;
        OX.state = otherX.state ;
        OX.invWW = otherX.invWW ;
        OX.N_cones=otherX.N_cones ;
        OX.ll    = otherX.ll ;
        OX.diff  = otherX.diff ;
        OX.beta  = otherX.beta ;
        OX.delta = otherX.delta ;
    else
        XX = X ;
        OX = otherX ;
    end
    Oclass   = classes{m}(classes{m}<=N) ;
    Xclass   = classes{m}(classes{m} >N) - N ;
    
%     fprintf('\t %d,%d',numel(Xclass),numel(Oclass)) ;
        
    if X.maxcones >= X.N_cones      + numel(Oclass)    && ...
       X.maxcones >= otherX.N_cones + numel(Xclass)

            xcx      = 1+mod(Xclass-1,X.M0) ;
            xcy      = 1+floor((Xclass-1)/X.M0) ;
            xcc      = XX.state(Xclass) ;
            
            ocx      = 1+mod(Oclass-1,X.M0) ;
            ocy      = 1+floor((Oclass-1)/X.M0) ;
            occ      = OX.state(Oclass) ;
   
            for k=1:length(Xclass)
                XX = flip_LL( XX , [xcx(k) xcy(k) 0] , PROB, T  ) ;
            end
            for k=1:length(Oclass)
                OX = flip_LL( OX , [ocx(k) ocy(k) 0] , PROB, oT ) ;
            end
            for k=1:length(Oclass)
                XX = flip_LL( XX , [ocx(k) ocy(k) occ(k)] , PROB, T  ) ;
            end
            for k=1:length(Xclass)
                OX = flip_LL( OX , [xcx(k) xcy(k) xcc(k)] , PROB, oT ) ;
            end
            
            XS{m}.X     = XX ;
            XS{m}.with  = OX ;
            XS{m}.ll    = XX.ll + OX.ll ;
            
        keep = [keep m] ;
    end
end
XS = XS(keep) ;


N = length(XS) ;
for i=1:N
    XS{i}.forward_prob = 1/N ;
end

if N==0
    X
    otherX
    keep
    numel(classes)
    classes{1}
    classes
    save(sprintf('crash_dump_%d.mat',randi(100)),'X','otherX','keep','classes','PROB')
end

end