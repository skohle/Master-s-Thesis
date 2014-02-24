function pkpod = adjointeqn_pod(n,deltat, yk, yqpod, c, ay, Mpod, Kpod, POD, Mpod4, M)

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
% Numerical methods:    POD-method,
%                       Semi-implicit Euler Method
%              
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Input parameters:
%
% y             % given POD state y
% m             % number of finite elements in [0,L]
% n             % number of time instances in [0,T]
% ydpod         % term yd of objective functional transf. with POD
% ay            % term ay of objective functional
% delta         % parameter delta in PDE   
% c             % parameter c in PDE
% Mpod,Kpod     % POD matrices
% POD           % POD basis
% Mpod4         % POD matrix phi^4
% deltat        % step size in time
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Output:   POD solution pkpod of the adjoint equation
% =========================================================================

[m,r] = size(POD);

[i,j,s] = find(POD'*M*ay);
ip = i;
jp = j+n-1;
sp = s;
ppost = sparse(i,j,s,r,1);



for k = n-1:-1:1
    Mck = zeros(r,r);
    for kp = 1:r
        for kq = 1:r
            Mck(kp,kq) = yk(:,k)'*Mpod4(:,:,kp,kq)*yk(:,k);
        end
    end
   
    % p(:,k) = (M/deltat + K + diag(exponent*c*M*(y(:,k).^(exponent-1))))\(M/deltat*p(:,k+1) + M*(y(:,k)-yqvec(:,k)));
	[i,j,s] = find((Mpod/deltat + Kpod + 3*c*Mck)\(Mpod/deltat*ppost + Mpod*(yk(:,k)-yqpod(:,k))));
    ip = [ip;i];
    jp = [jp;j+k-1];
    sp = [sp;s];
    ppost = sparse(i,j,s,r,1);
end
pkpod=sparse(ip,jp,sp,r,n);





