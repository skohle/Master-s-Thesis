function value=objvalue_pod(uk,yk,yqpod,ay,lambda,M,Mt,POD, Mpod)

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
% =====================

[r,n] = size(yqpod);

diffy = yk - yqpod;
value = 0.5*norm_Q(diffy, Mpod, Mt)^2 + (POD'*M*ay)'*Mpod*yk(:,n) + lambda/2*(uk(:,1)'*Mt*uk(:,1)+uk(:,2)'*Mt*uk(:,2)+uk(:,3)'*Mt*uk(:,3)+uk(:,4)'*Mt*uk(:,4));
