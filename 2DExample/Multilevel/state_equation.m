function y = state_equation(n,m,u,y0,c,fxyt,au,M,K,deltat)

% =========================================================================
% Author: Sophia Kohle, Technische Universit�t Berlin
% =========================================================================
%
% Task: Solution of the state equation to given controls u1,...,u4
%
%       y_t(x,t) - /\y(x,t) + c y^3(x,t) + d(x,t) = sum(i=1:4)au_i(x) u_i(t) 
%       y(x,0) = y_0(x) 
%       dy/dnu(x,t) = 0
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Numerical methods:    Finite Element Method,
%                       Semi-implicit Euler Method
%              
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Input parameters:
%
% u         % given control u
% y0        % initial value
% c         % parameter c in PDE
% fxyt      % parameter f(x,t) in PDE 
% au        % parameter au in PDE   
% M,K       % FE matrices
% deltat    % step size in time
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Output:   solution y of the state equation
% =========================================================================

[~,numcontrols] = size(u);
%m = length(y0);

% semi-implicit equation:
% M y'(t) + K y(t) + M f(x,t) + F(y(x,t)) = G * (u(t) + bt(t))
% y(0) = y0
[iy,jy,sy] = find(y0);
ypre = sparse(iy,jy,sy,m,1);
for k = 2:n
    % y(:,k) = (M/deltat + K)\((M*y(:,k-1))/deltat + sum(M*ai*ui(k)) - c*M*y(:,k-1).^3 - M*f(:,k));
    Muk = 0;
    for i = 1:numcontrols
        Muk = Muk + M*au(:,i)*u(k,i);
    end
    [i,j,s] = find((M/deltat + K)\((M*ypre)/deltat + Muk - c*M*ypre.^3 - M*fxyt(:,k)));
    iy = [iy,i];
    jy = [jy,j+k-1];
    sy = [sy,s];
    ypre = sparse(i,j,s,m,1);
end
y = sparse(iy,jy,sy,m,n);
