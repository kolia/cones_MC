function X = update_swap(trials,i)

if i>1

%     check_X(trials{1}.X)
%     check_X(trials{1}.with)

    X.X         = trials{1}.X ;
    X.X.state   = trials{i}.X.state ;
    if isfield(trials{i}.X,'invWW')
        X.X.invWW   = trials{i}.X.invWW ;
    end
    if isfield(trials{i}.X,'WW')
        X.X.WW      = trials{i}.X.WW ;
    end
    if isfield(trials{i}.X,'dUW_STA')
        X.X.dUW_STA = trials{i}.X.dUW_STA ;
    end
    if isfield(trials{i}.X,'ds_UW_STA')
        X.X.ds_UW_STA = trials{i}.X.ds_UW_STA ;
    end
    X.X.sparse_STA_W_state = trials{i}.X.sparse_STA_W_state ;
    X.X.N_cones = trials{i}.X.N_cones ;
    X.X.ll      = trials{i}.X.ll ;
    X.X.diff    = trials{i}.X.diff ;
    X.X.beta    = trials{i}.X.beta ;
    X.X.delta   = trials{i}.X.delta ;
    
    X.with         = trials{1}.with ;
    X.with.state   = trials{i}.with.state ;
    if isfield(trials{i}.with,'invWW')
        X.with.invWW   = trials{i}.with.invWW ;
    end
    if isfield(trials{i}.with,'WW')
        X.with.WW      = trials{i}.with.WW ;
    end
    if isfield(trials{i}.with,'dUW_STA')
        X.with.dUW_STA = trials{i}.with.dUW_STA ;
    end
    if isfield(trials{i}.with,'ds_UW_STA')
        X.with.ds_UW_STA = trials{i}.with.ds_UW_STA ;
    end
    X.with.sparse_STA_W_state = trials{i}.with.sparse_STA_W_state ;
    X.with.N_cones = trials{i}.with.N_cones ;
    X.with.ll      = trials{i}.with.ll ;
    X.with.diff    = trials{i}.with.diff ;
    X.with.beta    = trials{i}.with.beta ;
    X.with.delta   = trials{i}.with.delta ;
    
%     check_X(X.X)
%     check_X(X.with)
    
else
    X = trials{1} ;
end

% update both X
ii = 2*(i>1) + (i==1) ;
X.X     = update_X({trials{1}.X    X.X   },ii) ;
X.with  = update_X({trials{1}.with X.with},ii) ;

if isfield(X.X,'swap')
    X.X.swap(X.X.iteration) = true ;
end

if isfield(X.X,'swap')
    X.with.swap(X.with.iteration) = true ;
end

end