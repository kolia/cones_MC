function [classes,sizes] = equivalence_classes_direct(R,max_classes)
% for binary relation R, a boolean matrix
% calculate equivalence classes, i.e. groups of indices which are related
% through the transitive closure of R, as well as the sizes of these
% groups.


classes =  cell(max_classes,1) ;
sizes   = zeros(max_classes,1) ;

N_x = size(R,1) ;

inds_x  = find(sum(R,1)>0) ;
inds_y  = find(sum(R,2)>0) ;

R = R(inds_x,inds_y) ;

seen_x = zeros(numel(inds_x),1) ;
seen_y = zeros(numel(inds_y),1) ;

if ~isempty(inds_y)
    perm = randperm(numel(inds_y)) ;
    
    i=perm(1) ;
    k=1 ;
    while ~isempty(i)
        new_inds = i ;
        seen_y(i) = k ;
        while ~isempty(new_inds)
            s_x = sum(R(:,new_inds),2)>0 ;
            seen_x(s_x) = k ;
            s_y = sum(R(s_x,:),1)>0 ;
            new_inds = find(s_y' .* (seen_y>0)) ;
            seen_y(new_inds) = k ;
        end

        s_x = seen_x==k ;
        s_y = seen_y==k ;
        
        classes{k} = [inds_x(s_x) ; N_x+inds_y(s_y)] ;
        sizes(k)   = numel(classes{k}) ;

        if k >= max_classes
            break
        end
        k = k+1 ;
        i = perm( find(~seen_y(perm),1) ) ;
    end
    classes = classes(1:k) ;
    sizes   = sizes(1:k) ;
else
    classes = {} ;
    sizes   = [] ;
end
end