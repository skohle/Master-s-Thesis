function p = adjoint_equation(n,m,y,deltat,yq,c,ay,M,K)

% =========================================================================
% Author: Sophia Kohle, Technische Universität Berlin
% =========================================================================
%
% Task: Solution of the adjoint equation to a given state y
%
%       -p_t(x,t) - /\p(x,t)+3c*y^2(x,t)p(x,t)+d(x,t) = y(x,t) - yq(x,t)   in Q                               in (0,L) x (0,T)
%       p(x,T) = ay;                                   in Omega
%       dp/dnu(x,t)=0           in Sigma
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
% y         % given state y
% m         % number of finite elements in [0,L]
% n         % number of time instances in [0,T]
% yd        % term yd of objective functional
% ay        % term ay of objective functional
% delta     % parameter delta in PDE   
% c         % parameter c in PDE
% M,K,Q,G   % FE matrices
% deltat    % step size in time
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Output:   solution p of the adjoint equation
% =========================================================================

% implicit equation:
% -M p'(t) + K p(t) + delta Q y(t) + 3c G y(x,t)^2p(x,t) = -ay G
% p(T) = ay
[i,j,s] = find(ay);
ip = i;
jp = j+n-1;
sp = s;
ppost = sparse(i,j,s,m,1);
for k = n-1:-1:1
    
    % p(:,k) = (M/deltat + K + diag(exponent*c*M*(y(:,k).^(exponent-1))))\(M/deltat*p(:,k+1) + M*(y(:,k)-yqvec(:,k)));
	[i,j,s] = find((M/deltat + K + diag(3*c*M*(y(:,k).^2)))\(M/deltat*ppost + M*(y(:,k)-yq(:,k))));
    ip = [ip;i];
    jp = [jp;j+k-1];
    sp = [sp;s];
    ppost = sparse(i,j,s,m,1);
end
p = sparse(ip,jp,sp,m,n);
