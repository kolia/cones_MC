function evidence = plotable_evidence( evidence )

z = sum(evidence,3) ;
inds  = logical(z<0) ;
inds  = inds(:) ;
inds3 = [inds ; inds ; inds] ;
evidence(inds3) = -[z(inds) ; z(inds) ; z(inds)] ;

evidence = evidence-min(evidence(:)) ;
evidence = evidence / max(evidence(:)) ;

end