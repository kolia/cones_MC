function [x,y,states] = get_best( Xs )
% Xs: cell struct with either bestX or X field.
% Each bestX (or X) is either a struct, in which case
% x and y are taken from the 'N_cones' and 'll' fields resp.
% or bestX (or X) is a cell, in which x and y are taken from the same
% fields of the bestX (or X) with the largest 'll'.

x      = zeros(numel(Xs),1) ;
y      = zeros(numel(Xs),1) ;
states = cell(numel(Xs),1) ;
for i=1:numel(Xs)
    if isfield(Xs{i},'bestX')
        X = Xs{i}.bestX ;
    else
        X = Xs{i}.X ;
    end
    if isa(X, 'cell')
        j = 0 ;
        bestll = -1e10 ;
        for k=1:numel( X )
            if X{k}.ll > bestll
                j = k ;
            end
        end
        x(i) = X{j}.N_cones ;
        y(i) = X{j}.ll ;
        states{i} = X{j}.state ;
    else
        x(i) = X.N_cones ;
        y(i) = X.ll ;
        states{i} = X.state ;
    end
end

end