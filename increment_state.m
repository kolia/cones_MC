function state = increment_state( dX , state )

n = numel(max(dX(:,1:3:end))) ;
N = find(max(dX,[],2),1,'last') ;

for i=1:N
    for j=0:n-1
        if dX(i,1+3*j)
            state(dX(i,1+3*j),dX(i,2+3*j)) = dX(i,3+3*j) ;
        end
    end
end

end