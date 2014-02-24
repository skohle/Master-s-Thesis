function pkpod = adjointeqn_pod_deim(n,deltat, yk, yqpod, c, ay, Mpod, Kpod, POD, M,U,P)

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
% yk            % given POD state y
% m             % number of finite elements in [0,L]
% n             % number of time instances in [0,T]
% ydpod         % term yd of objective functional
% ay            % term ay of objective functional
% delta         % parameter delta in PDE   
% c             % parameter c in PDE
% Mpod,Kpod     % POD matrices
% POD           % POD basis
% M             % massmatrix FEM
% U,P           % POD-DEIm matrices
% deltat        % step size in time
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Output:   POD solution p of the adjoint equation
% =========================================================================

[m,r] = size(POD);

[i,j,s] = find(POD'*M*ay);
ip = i;
jp = j+n-1;
sp = s;
ppost = sparse(i,j,s,r,1);

Proj=(POD'*M*U)/(P'*U);
PPOD=P'*POD;

for k = n-1:-1:1
    
    % p(:,k) = (M/deltat + K + diag(exponent*c*M*(y(:,k).^(exponent-1))))\(M/deltat*p(:,k+1) + M*(y(:,k)-yqvec(:,k)));
	[i,j,s] = find((Mpod/deltat + Kpod + diag(3*c*Proj*(PPOD*yk(:,k)).^2))\(Mpod/deltat*ppost + Mpod*(yk(:,k)-yqpod(:,k))));
    ip = [ip;i];
    jp = [jp;j+k-1];
    sp = [sp;s];
    ppost = sparse(i,j,s,r,1);
end
pkpod=sparse(ip,jp,sp,r,n);






