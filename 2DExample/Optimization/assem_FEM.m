function [M,K,F,Q,G,H,R] = assem_FEM(xy,edges,tri,gammashape,au)

% =========================================================================
% Author: Sophia Kohle, Technische Universität Berlin
% =========================================================================
%
% Task: Matrix assembling (linear FE ansatz functions)
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Input parameters:
%              
% xy,edges,tri      % nodes, edges, triangles of FEM in space
% gammashape        % boundary file 
% au                % term au of objective functional
%                
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Output:   matrices M, K, Q, G, F, H, R
% =========================================================================

for i = 1:size(au,2)
    autri(:,i) = pdeintrp(xy,tri,au(:,i));
    [K,M,F(:,i)] = assema(xy,tri,1,1,autri(:,i)');
end

[Q,G,H,R] = assemb(gammashape,xy,edges);
K = sparse(K);
M = sparse(M);
F = sparse(F);
Q = sparse(Q);
G = sparse(G);
H = sparse(H);
R = sparse(R);
