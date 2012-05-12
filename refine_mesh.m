function [Y,Z] = refine_mesh(X, curvature)
% Consider each column of X to be the values of a function on a mesh.
% Refine this mesh by making X larger in the first dimension, adding values
% that interpolate linearly between the old ones.  The first dimension of
% the result is 2 * the old one - 1.
% If X is a cell, then every element in the cell is treated as a row of X,
% and the result is 
%

% if X is a cell, convert to matrix
Xcell = iscell(X) ;
if Xcell, X = cell2mat(X(:)) ; end

Y = zeros(2*size(X,1)-1,size(X,2)) ;
Y(1:2:end,:) = X ;
Y(2:2:end,:) = X(1:end-1,:).*(1-curvature)+curvature.*X(2:end,:) ;

Z = Y ;

% if X was a cell, convert result back to cell
if Xcell
    Y = cell(1,size(X,1)) ;
    for i=1:size(Z,1)
        Y{i} = Z(i,:) ;
    end
end
end