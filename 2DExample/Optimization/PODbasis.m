function [ysnap,POD] = PODbasis(n,m,r,deltat,c,fxyt,au,usnap,M,K,Mt,y0,ysnap,mode_d)

% =========================================================================
% Author: Sophia Kohle, Technische Universität Berlin
% =========================================================================
%
% Task: Computation of POD basis of rank r for the state obtained as
%       solution of
%
%       y_t(x,t) - /\y(x,t) + c y^3(x,t) + d(x,t) = sum(i=1:4)au_i(x) u_i(t) 
%       y(x,0) = y_0(x) 
%       dy/dnu(x,t) = 0
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Input parameters:
%
% r         % rank of POD basis
% deltat    % step size in time
% c         % parameter c in PDE
% fxyt      % term fxyt on left-hand side of PDE 
% au        % term au of objective functional
% usnap     % given control to generate snapshots
% M,K       % FE matrices
% Mt        % time discr. matrix
% y0        % initial value
% ysnap     % snapshots (case POD-DEIM)
% mode_d    % 0: snapshots for whole PDE
%           % 1: snapshots for nonlinearity
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Output:   POD basis, Snapshots
% =========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% snapshots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mode_d==0 %for the whole PDE
    clear ysnap
    
    ysnap = state_equation(n,m,usnap,y0,c,fxyt,au,M,K,deltat);
    disp('--- Start Computing POD basis ------------------------------------');
else %for the nonlinearity if compitation with POD-DEIM
    ysnap=ysnap.^3;
    disp('--- Start Computing POD basis for the nonlinear term -------------');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% singular value decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mtsqrt = sqrt(Mt);
YMY = Mtsqrt*ysnap'*M*ysnap*Mtsqrt;
YMY = 0.5*(YMY + YMY'); % for symmetry reasons
[EV,EW] = eig(YMY);

% choose r biggest eigenvalues
DEW = diag(EW);
maxEW = sparse(r,1);
indmax = sparse(r,1);
for i = 1:r
    [maxEW(i),indmax(i)] = max(DEW);
    DEW(indmax(i)) = -Inf;
end

    disp('Decay of eigenvalues in POD basis computation:')
    disp(' ');
    for i = 1:r
        disp(['   ' num2str(maxEW(i))]);
    end

% eigenvectors
V = EV(:,indmax);

% format
oversqrtmaxEW = (1./(sqrt(maxEW)));
VoversqrtEW = V * diag(oversqrtmaxEW);

% POD basis vector transformation
POD = ysnap*Mtsqrt*VoversqrtEW;