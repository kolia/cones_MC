function [results,X] = update_X(results,trials,i)

changed             = i>1 ;
X                   = trials{i} ;
if i>1 , X.version  = trials{1}.version+1 ; end
results.iteration   = results.iteration+1 ;
clear trials


%% accumulate statistics

if ~isfield(X,'burn_in')  % only after burn_in iterations    
    
    % accumulate summed state
    if isfield(results,'summed')
        if changed
            results.summed = results.summed + results.change * results.times ;
            results.change = ...
                logical( [(X.state(:)'==1) (X.state(:)'==2) (X.state(:)'==3)] ) ;
            results.times = 1 ;
        else
            results.times = results.times + 1 ;
        end
        results.N_iter  = results.N_iter + 1 ;
    end
    
    if isfield(results,'LL_mean')
        results.N500     = (1-1/500)*results.N500     + 1 ;
        results.LL_mean  = (1-1/500)*results.LL_mean  + X.ll ;
    end
    
    if isfield(results,'LL_mean_square')
        results.LL_mean_square = (1-1/500)*results.LL_mean_square + X.ll^2 ;
    end
    
end

%% track all changes
if isfield(results,'dX')
    if changed
        diffX = X.diff' ;
        results.dX(results.iteration,1:numel(diffX)) = diffX(:)' ;
    else
        results.dX(results.iteration,1) = 0 ;
    end
end


if changed
        
    %% update contacts

    deleted = find( ~X.diff(:,3) ) ;
    if ~isempty(deleted)
        deleted = X.diff(deleted,1) + (X.diff(deleted,2)-1)*X.M0 ;
        
        for d=1:4
            X.contact{d}(deleted,:) = false ;
            X.contact{d}(:,deleted) = false ;
        end
    end
    
    added   = find( X.diff(:,3) ) ;
    for d=1:4
        for i=1:numel(added)
            X = make_contacts(X , X.diff(added(i),1) , X.diff(added(i),2)) ;
        end
    end
    
%     check_X(X)
    
    X.diff    = [] ;
    
end
end