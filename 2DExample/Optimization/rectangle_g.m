function [geometry] = rectangle_g(x1,x2,y1,y2)

% =========================================================================
% Author: Eileen Kammann, Technische Universität Berlin
% =========================================================================
%
% Geometry File defining the geometry of a rectangle. 
% Used for initmesh
%
% =========================================================================

geometry = [2 2 2 2             % 2 = line segment; 1 = circle segment; 4 = elipse segment
            x1 x2 x2 x1         % start point x
            x2 x2 x1 x1         % end point x
            y1 y1 y2 y2         % start point y
            y1 y2 y2 y1         % end point y
            1 1 1 1 
            0 0 0 0];