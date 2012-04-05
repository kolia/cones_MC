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

if ~isempty(X.diff)
        
    %% update contacts
        
    deleted = find( ~X.diff(:,3) ) ;
    if ~isempty(deleted)
        if track_contacts && isfield(X,'contact')
            for ii=1:numel(deleted)
                [~,indices] = not_excluded( X, X.diff(deleted(ii),1), X.diff(deleted(ii),2) ) ;
                X.excluded(indices) = false ;
            end

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
    end

    added   = find( X.diff(:,3) ) ;
    for ii=1:numel(added)
        if track_contacts && isfield(X,'contact')
            X = make_contacts(X , X.diff(added(ii),:)) ;
        end
        [~,indices] = not_excluded( X, X.diff(added(ii),1), X.diff(added(ii),2) ) ;
        X.excluded(indices) = true ;
    end
end


%% track all changes
if isfield(X,'LL_history')
    if X.iteration>length(X.LL_history)
       X.LL_history = [X.LL_history ; zeros(numel(X.LL_history)+500,1)] ;
    end
    X.LL_history(X.iteration) = X.ll ;
end

if isfield(X,'N_cones_history')
    if X.iteration>length(X.N_cones_history)
       X.N_cones_history = [X.N_cones_history ; zeros(numel(X.N_cones_history)+500,1)] ;
    end
    X.N_cones_history(X.iteration) = X.N_cones ;
end

if X.iteration>numel(X.cputime)
    X.cputime = [X.cputime ; zeros(numel(X.cputime)+500,1)] ;
end
X.cputime(X.iteration) = cputime ;
X.diff    = [] ;
    
end