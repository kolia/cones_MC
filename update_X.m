function X = update_X(trials,i,track_contacts)

if nargin<2 ,  i = 1 ; end
if nargin<3 ,  track_contacts=true ; end
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

%% track all changes
if isfield(X,'LL_history')
    if X.iteration>length(X.LL_history)
       X.LL_history = [X.LL_history ; zeros(500,1)] ;
    end
    X.LL_history(X.iteration) = X.ll ;
end

if track_contacts && ~isempty(X.diff)
        
    %% update contacts
    if isfield(X,'contact')
        
        deleted = find( ~X.diff(:,3) ) ;
        if ~isempty(deleted)
            deleted = X.diff(deleted,1) + (X.diff(deleted,2)-1)*X.M0 ;
            
            for d=1:2
                X.contact{d}(deleted,:) = false ;
%                 test = X.contact{d} ;
%                 test(deleted,:) = false ;
%                 for del=deleted
%                     inds = X.contact{d}(del,:) ;
%                     X.contact{d}(del,inds) = false ;
%                 end
%                 try X.contact{d} - test ;
%                 catch                    
%                     'oups'
%                 end
                X.contact{d}(:,deleted) = false ;
            end
        end
        
        added   = find( X.diff(:,3) ) ;
        for i=1:numel(added)
            X = make_contacts(X , X.diff(added(i),:)) ;
        end
    end
end

if X.iteration>numel(X.cputime)
    X.cputime = [X.cputime ; zeros(numel(X.cputime,1))] ;
end
X.cputime(X.iteration) = cputime ;
X.diff    = [] ;
    
end