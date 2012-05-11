function [slowX,fastX] = update_swap(slowX,fastX,proposed_slowX,proposed_fastX,accept)

% update both X
slowX  = update_X({slowX proposed_slowX},accept+1) ;
fastX  = update_X({fastX proposed_fastX},accept+1) ;

if isfield(slowX,'swap')
    slowX.swap(slowX.iteration) = true ;
end

if isfield(fastX,'swap')
    fastX.swap(fastX.iteration) = true ;
end

end