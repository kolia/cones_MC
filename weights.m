function [m,C,scaling,X] = weights( x , PROB )
% Get mean m and covariance of weights between cones and ganglion cells
% Weights between a cone and different ganglion cells are independent.
% The covariance matrix between cone weights is the same for every ganglion
% cell, up to a scaling factor: the covariance between cones i and j for
% ganglion cell g is given by  C(i,j)*scaling(g).
% 
% INPUT:
%   x           either a 2D matrix with 0,1,2 or 3 in each entry
%               where 0 is no cone, 1 = red, 2=green, 3=blue
%               or a set of incremental moves, specified by an N-by-3
%               matrix of the form [x y c], where the x column contains x
%               coordinates of cones, the y column contains y coordinates,
%               and c contains colors (1,2 or 3).
%
%   PROB        struct containing problem specification.
%               returned by exact_LL_setup.m and cones.m
%               required fields:
%               - M0
%               - M1 
%               - SS
%               - coneConv
%               - colorDot
%               - STA_W
%               - sumLconst
%               - maxcones
%               - N_colors
%               - D

M0 = cone_map.M0 * cone_map.cone_params.supersample ;
M1 = cone_map.M1 * cone_map.cone_params.supersample ;

%  if an X structure is given instead of a x,
%  extract the x first
if isstruct(x) && isfield(x,'state')
    x = x.state ;
end

%  if a set of moves was given instead of a x
%  then calculate x from moves
if max(x(:)) > 3
    x=increment_state(x,zeros(M0,M1) ) ;
end

[x,y,c] = find(x) ;

dX = [x y c] ;

X = initialize_X( PROB.M0, PROB.M1, PROB.N_colors, PROB.SS, PROB.D, PROB.naive_LL, 1 ) ;

X = flip_LL( X , dX , PROB ) ;
    
STA_W_state = PROB.STA_W( x+M0*(y-1)+M0*M1*(c-1) , : ) ;


m = X.invWW * STA_W_state ;
C = inv( X.invWW ) ;
scaling = 1 ./ PROB.cell_consts ;

end