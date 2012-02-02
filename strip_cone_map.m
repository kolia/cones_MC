for i=1:50
    try
        name   = sprintf('result_%d.mat',i) ;
        load( name ) ;
        try
            cone_map = rmfield(cone_map,{'STA','LL','NICE','initX'}) ;
        end
        cone_map.X = rmfield(cone_map.X,{'invWW','contact'}) ;
        save(name,'cone_map')
    end
end

names = fieldnames( cone_map ) ;                                          
for i=1:numel(names) , assignin('base',names{i},cone_map.(names{i})) ; end