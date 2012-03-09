function [classes,sizes] = equivalence_classes(R,max_classes)
% for binary transitive relation R, a boolean square matrix
% calculate equivalence classes, i.e. groups of indices which are related,
% as well as the sizes of these groups.


n    = size(R,1) ;
seen = zeros(n,1) ;

classes =  cell(min(n,max_classes),1) ;
sizes   = zeros(min(n,max_classes),1) ;

inds    = find(diag(R)>0) ;

if ~isempty(inds)
    perm = randperm(numel(inds)) ;
    inds = inds(perm) ;
    
    i=inds(1) ;
    k=1 ;
    while ~isempty(i)
        s = R(:,i) ;
        seen(s) = k ;

        classes{k} = find(s) ;
        sizes(k)   = sum(s) ;

        if k >= max_classes
            break
        end
        k = k+1 ;
        i = inds( find(~seen(inds),1) ) ;
    end
    classes = classes(1:k) ;
    sizes   = sizes(1:k) ;
else
    classes = {} ;
    sizes   = [] ;
end
end