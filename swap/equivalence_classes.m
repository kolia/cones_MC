function [classes,sizes] = equivalence_classes(R)
% for binary transitive relation R, a boolean square matrix
% calculate equivalence classes, i.e. groups of indices which are related,
% as well as the sizes of these groups.


n    = size(R,1) ;
seen = zeros(n,1) ;

classes =  cell(n,1) ;
sizes   = zeros(n,1) ;

inds    = find(diag(R)>0) ;

i=inds(1) ;
k=1 ;
while ~isempty(i)
    s = R(i,:) ;
    seen(s) = k ;
    
    classes{k} = find(s) ;
    sizes(k)   = sum(s) ;
    
    k = k+1 ;
    
    inds = inds(inds>i) ;
    
    i = inds( find(~seen(inds),1) ) ;
end

classes = classes(1:k-1) ;
sizes   =   sizes(1:k-1) ;

end