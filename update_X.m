function [results,X] = update_X(results,X,changed)


%% accumulate statistics

if ~isfield(X,'burn_in')  % only after burn_in iterations    
    
    % accumulate summed state
    if isfield(results,'summed')
        if changed
            results.summed = results.summed + results.change * results.times ;
            results.change = ...
                [(X.state(:)'==1) (X.state(:)'==2) (X.state(:)'==3)] ;
            results.times = 1 ;
        else
            results.times = results.times + 1 ;
        end
        results.N_iter  = results.N_iter + 1 ;
    end
    

end

if changed
    %% update contacts
    
    
    if ~isempty(X.diff.deleted)
        deleted = X.diff.deleted(:,1) + (X.diff.deleted(:,2)-1)*X.M0 ;
        
        for d=1:4
            X.contact{d}(deleted,:) = false ;
            X.contact{d}(:,deleted) = false ;
        end
    end
    
    
    for d=1:4
        for i=1:size(X.diff.added,1)
            X = make_contacts(X , X.diff.added(i,1) , X.diff.added(i,2)) ;
        end
    end
    
    % check_X(X)
    
    X.diff.added    = [] ;
    X.diff.deleted  = [] ;
    
end
end