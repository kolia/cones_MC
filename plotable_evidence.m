function evidence = plotable_evidence( evidence )

infinite = ~isfinite(evidence) ;
evidence(~isfinite(evidence)) = 1e-6 ;

z = max(evidence,[],3) ;
inds  = logical(z<=0) ;
inds  = inds(:) ;
inds3 = [inds ; inds ; inds] ;
% evidence(inds3) = -[z(inds) ; z(inds) ; z(inds)] ;
evidence(inds3) = 0 ;

evidence = evidence-min(evidence(:)) ;
evidence = evidence / max(evidence(:)) ;

evidence(inds3) = 0.001 ;
evidence(infinite) = 0 ;

end