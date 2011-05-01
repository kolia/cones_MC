function X = update_X(trials,i)

if nargin<2 ,  i = 1 ; end
changed             = i>1 ;
X                   = trials{i} ;
if i>1 , X.version  = trials{1}.version+1 ; end
if ~isfield(X,'iteration') , X.iteration = 0 ; end
X.iteration   = X.iteration+1 ;
clear trials


%% track all changes
if isfield(X,'dX')
    if changed
        diffX = X.diff' ;
        X.dX(X.iteration,1:numel(diffX)) = diffX(:)' ;
    else
        X.dX(X.iteration,1) = 0 ;
    end
end


if changed
        
    %% update contacts
    if isfield(X,'contact')
    
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
                X = make_contacts(X , X.diff(added(i),:)) ;
            end
        end
        
%             check_X(X)        
    end    
end

X.diff    = [] ;
    
end