function [boundary] = rectangle_b(stringb1, stringb2, stringb3, stringb4)

% =========================================================================
% Author: Eileen Kammann, Technische Universität Berlin
% =========================================================================
%
% Boundary File defining the boundary conditions of a rectangle. 
% Used for assemb
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% du/dny (x) + q u = g
%              h u = r
%
%   1: #PDEs
%   2: #Dirichlet boundary conditions
%   3: length q
%   4: length g
%   7: sting q
%   8: string g
%  (5: length h)
%  (6: length r)
%  (9: string h)
% (10: string r)
%
% =========================================================================

b1 = [1, 0, 1, length(stringb1), '0', stringb1]';
b2 = [1, 0, 1, length(stringb2), '0', stringb2]';
b3 = [1, 0, 1, length(stringb3), '0', stringb3]';
b4 = [1, 0, 1, length(stringb4), '0', stringb4]';

boundary = [ b1, b2, b3, b4 ];