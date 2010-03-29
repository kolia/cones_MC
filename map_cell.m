function data = map_cell( f , data )

if iscell(data)
    [N,M] = size(data) ;    
    for i=1:N
        for j=1:M
            data{i,j} = f( data{i,j} ) ;
        end
    end
end

end