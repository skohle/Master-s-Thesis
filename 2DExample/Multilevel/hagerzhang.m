function [betak]=hagerzhang(gk,gl,dk)
% =========================================================================
% Author: Sophia Kohle, Technische Universität Berlin
% =========================================================================
% 
% Update-formula Hager-Zhang
%
% Input:
%
% gk    % current gradient
% gl    % gradient step k-1
% dk    % decent
%
% Output:
%
% betak     Update
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta=0.01;
gamma=gk-gl;
beta=(gamma-2*(gamma'*gamma)^2/(dk'*gamma)*dk)'*(gk/(dk'*gamma));
etak=-1/(sqrt(dk'*dk))*min(eta,sqrt(gl'*gl));
betak=max(beta,etak);
