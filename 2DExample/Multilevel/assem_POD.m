function [Mpod,Kpod,y0pod,yqpod] = assem_POD(M,K,y0,yq,POD)
         
% =========================================================================
% Author: Sophia Kohle, Technische Universität Berlin
% =========================================================================
%
% Task: Matrix assembling (given POD basis functions = ansatz functions)
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Input parameters:
%              
% m         % number of finite elements in [0,L]
% yd        % term yd of objective functional
% y0        % initial value
% M,K       % FE matrices
% POD       % POD basis vectors
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Output:   POD matrices Mpod, Qpod, Kpod, Gpod, y0pod, ydpod
% =========================================================================

r = size(POD,2);

% POD matrix assembling
% ------------------------------------------------------------------------
% Mpod = POD'*M*POD;
Mpod = speye(r);

Kpod = POD'*K*POD;
Kpod = 0.5.*(Kpod + Kpod'); % eliminate rounding error

y0pod = Mpod*POD'*M*y0;

yqpod = POD'*M*yq;


