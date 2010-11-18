function [init_X , dX] = load_stats()

lsd = ls('stats_*.mat') ;
filenames = {} ;

if ~strcmp(lsd,'ls: ')  % files exist
    filenames = regexp(lsd,'\s','split') ;
end

dX = cell( length(filenames) , 1 ) ;
init_X = cell( length(filenames) , 1 ) ;
ii = 1 ;
for i=1:length(dX)
    if ~isempty(filenames{i})
        s = load(filenames{i}) ;
        dX{ii} = full( s.results{1}.dX ) ;
        init_X{ii} = s.initial_X{1}.state ;
        ii = ii + 1 ;
    end
end
init_X = init_X(1:ii-1) ;
dX = dX(1:ii-1) ;

end