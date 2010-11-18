function cone_map = coalescer( cone_map , folder , pattern , instances , parameters , ID )

here = pwd ;
cd(folder)

filenames = regexp( ls , pattern , 'match' ) ;

cone_map.X = {} ;

for i=1:length(instances)    
    bestX = load( filenames{instances(i)} ) ;
    
    if isstruct(bestX) && isfield(bestX,'bestX')
        bestX = bestX.bestX ;
    end
    
    if iscell(bestX) && length(bestX) == 1
        bestX = bestX{1} ;
    end
    
    if ~isfield( bestX,'ll')
        cone_map.X = [cone_map.X bestX.bestX] ;
    else
        cone_map.X = [cone_map.X bestX] ;
    end
end

cd(here)

if isfield(cone_map,'ROI')
    cone_map    = rmfield(cone_map,'ROI') ;
end

cone_map.ID = ID ;

cone_map = parameters( cone_map ) ;

end