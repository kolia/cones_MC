function XS = swap_closure( X , T, otherX , oT, PROB )
% Combine two systems into a single system, summing their log-likelihoods.
% This is in preparation for swapping parts of the two systems'
% configurations. For this, groups of cells that must be swapped together
% are calculated.

% check_X(X)
% check_X(otherX)

% calculate overlaps of X cones on otherX exclusion disks

% symmetrize relation and calculate transitive closure
N = numel(X.state) ;

% i = find( otherX.state) ;
% j = find( X.state) ;
R = overlap_relation( otherX , X ) ;

% EO= logical(sparse(i,i,ones(length(i),1),N,N,3*N)) ;
% EX= logical(sparse(j,j,ones(length(j),1),N,N,3*N)) ;
% R = logical( [ EO R ; R' EX ] ) ;
% R = transitive_closure(R,[i ; N+j],0) ;

% EO= logical(sparse([],[],[],N,N,N)) ;
% EX= logical(sparse([],[],[],N,N,N)) ;
% RR = logical( [ EO R ; R' EX ] ) ;
% RR = transitive_closure(RR,[i ; N+j],1) ;

% clear i j

% get equivalence classes / connected components
% classes_old = equivalence_classes(RR,20) ;

classes  = equivalence_classes_direct(R,50) ;

if numel(classes)>0
    classes = [{[]} ; classes] ;
else
    classes = {[]} ;
end
XS = cell( numel(classes), 1 ) ;

fprintf(' %d',numel(classes))


keep = [] ;
for m=1:length(XS)
    if m>1
        XX       = struct ;
        XX.state = X.state ;
        if isfield(X,'invWW')
            XX.invWW = X.invWW ;
        end
        if isfield(X,'WW')
            XX.WW = X.WW ;
        end
        if isfield(X,'dUW_STA')
            XX.dUW_STA = X.dUW_STA ;
        end
        if isfield(X,'ds_UW_STA')
            XX.ds_UW_STA = X.ds_UW_STA ;
        end
        XX.sparse_STA_W_state = X.sparse_STA_W_state ;
        XX.contributions = X.contributions ;
        XX.N_cones=X.N_cones ;
        XX.ll    = X.ll ;
        XX.diff  = X.diff ;
        XX.beta  = X.beta ;
        XX.delta = X.delta ;
        XX.LL_history = X.LL_history ;
        XX.N_cones_history = X.N_cones_history ;
        XX.cputime    = X.cputime ;
        
        OX       = struct ;
        OX.state = otherX.state ;
        if isfield(X,'invWW')
            OX.invWW = otherX.invWW ;
        end
        if isfield(X,'WW')
            OX.WW = otherX.WW ;
        end
        if isfield(X,'dUW_STA')
            OX.dUW_STA = otherX.dUW_STA ;
        end
        if isfield(X,'ds_UW_STA')
            OX.ds_UW_STA = otherX.ds_UW_STA ;
        end
        OX.sparse_STA_W_state = otherX.sparse_STA_W_state ;
        OX.contributions = otherX.contributions ;
        OX.N_cones=otherX.N_cones ;
        OX.ll    = otherX.ll ;
        OX.diff  = otherX.diff ;
        OX.beta  = otherX.beta ;
        OX.delta = otherX.delta ;
        OX.LL_history = otherX.LL_history ;
        OX.N_cones_history = otherX.N_cones_history ;
        OX.cputime    = otherX.cputime ;
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


N = numel(XS) ;
for i=1:N
%     check_X(XS{i}.X)
%     check_X(XS{i}.with)
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