function [m,xy,edges,tri] = grid_generation(mmin, omegashape)

% =========================================================================
% Author: Sophia Kohle, Technische Universität Berlin
% =========================================================================
%
% Discretization of some space Omega with geometry file omegashape by 
% Finite Elements Method
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Input parameters:
%              
% mmin       % minimal number of finite elements in space discretization
% omegashape % geometry file for space Omega
%                
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Output parameters:
% xy,edges,tri      % nodes, edges, triangles of FEM in space
%
% =========================================================================
% =========================================================================

[xy,edges,tri] = initmesh(omegashape,'Hgrad',1.05);
while size(xy,2) < mmin
    [xy,edges,tri] = refinemesh(omegashape,xy,edges,tri);
end
m = size(xy,2);

