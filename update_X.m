function X = update_X(trials,i,track_contacts)
% X = update_X(trials,[i,track_contacts])
% choose trials{i} as next X
% update X based on changes since last iteration recorded in X.diff
% update X.contact if track_contacts is true
% record X.LL_history, X.N_cones_history, X.cputime if they exist

% by default,  i = 1   which means the empty move was chosen
if nargin<2 ,  i = 1 ; end

% default is to track contacts
if nargin<3 ,  track_contacts=true ; end

% choose trials{i} as next X
X                   = trials{i} ;
if i>1 , X.version  = trials{1}.version+1 ; end

% update iteration number
if ~isfield(X,'iteration') , X.iteration = 0 ; end
X.iteration   = X.iteration+1 ;

% all other proposed moves can be discarded
clear trials


%% track all changes if dX exists (used to replay CAST runs later)
if isfield(X,'dX')
    if i>1
        diffX = X.diff' ;
        X.dX(X.iteration,1:numel(diffX)) = diffX(:)' ;
    else
        X.dX(X.iteration,1) = 0 ;
    end
end

%% update data structures used by greedy, MCMC and CAST
if ~isempty(X.diff)
    deleted = find( ~X.diff(:,3) ) ;
    if ~isempty(deleted)
        %% update excluded (used by all algorithms)
        for ii=1:numel(deleted)
            [~,indices] = not_excluded( X, X.diff(deleted(ii),1), X.diff(deleted(ii),2) ) ;
            X.excluded(indices) = false ;
        end
        
        %% update contacts (used by MCMC and CAST for shift moves)
        if track_contacts && isfield(X,'contact')

            deleted = X.diff(deleted,1) + (X.diff(deleted,2)-1)*X.M0 ;            
            for d=1:2
                X.contact{d}(deleted,:) = false ;
                X.contact{d}(:,deleted) = false ;
            end
        end
    end

    added   = find( X.diff(:,3) ) ;
    for ii=1:numel(added)
        %% update contacts (used by MCMC and CAST for shift moves)
        if track_contacts && isfield(X,'contact')
            X = make_contacts(X , X.diff(added(ii),:)) ;
        end
        %% update excluded (used by all algorithms)
        [~,indices] = not_excluded( X, X.diff(added(ii),1), X.diff(added(ii),2) ) ;
        X.excluded(indices) = true ;
    end
end


%% record LL_history if field exists
if isfield(X,'LL_history')
    if X.iteration>length(X.LL_history)
       X.LL_history = [X.LL_history ; zeros(numel(X.LL_history)+500,1)] ;
    end
    X.LL_history(X.iteration) = X.ll ;
end

%% record N_cones_history if field exists
if isfield(X,'N_cones_history')
    if X.iteration>length(X.N_cones_history)
       X.N_cones_history = [X.N_cones_history ; zeros(numel(X.N_cones_history)+500,1)] ;
    end
    X.N_cones_history(X.iteration) = X.N_cones ;
end

%% record cputime if field exists
if isfield(X,'cputime')
    if X.iteration>numel(X.cputime)
        X.cputime = [X.cputime ; zeros(numel(X.cputime)+500,1)] ;
    end
    X.cputime(X.iteration) = cputime ;
end

%% reset changes since last update to empty change set
X.diff    = [] ;
    
end