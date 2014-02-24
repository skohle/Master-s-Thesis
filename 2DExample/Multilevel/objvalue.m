function value = objvalue(ui,y,yq,ay,lambda,M,Mt)

% =========================================================================
% Author: Sophia Kohle, Technische Universität Berlin
% =========================================================================
%
% Task: Computation of the function value of 
%
%       J(y,ui) := 1/2 ||y(x,t)-y_Q(x,t)||^2_L^2(Q) + (ay(x), y(l,T))
%                 + lambda/2 sumi||ui(t)||^2_L^2(0,T)
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Input parameters:
%
% ui        % given controls ui
% y         % given state y
% yq        % parameter yq objective function
% ay        % parameter c in objective function
% lambda    % parameter lambda in objective function
% M         % FE mass matrix
% Mt
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Output:   J(y,ui)
% =========================================================================

[m,n] = size(yq);

diffy = y - yq;
value = 0.5*norm_Q(diffy, M, Mt)^2 + ay'*M*y(:,n) + lambda/2*(ui(:,1)'*Mt*ui(:,1)+ui(:,2)'*Mt*ui(:,2)+ui(:,3)'*Mt*ui(:,3)+ui(:,4)'*Mt*ui(:,4));
