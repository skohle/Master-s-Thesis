function [sk,uk,ykpod_deim,pkpod_deim,J_uk,gk] = linesearch_WP_pod_deim(r,ul,gl,J_ul,s0,dk,y0pod,c,fxyt,au,M,deltat,yqpod,ay,lambda,Mt,ua,ub,Mpod,Kpod,POD,V1pod,V2pod,V3pod,V4pod,U,P, x1x2,r_deim)
% =========================================================================
% Author: Sophia Kohle, Technische Universität Berlin
% =========================================================================
%
% Task: Computation of an optimal control for the problem
%
% min 1/2 ||y-yq||^2 + lambda/2 sum(i=1:4) ||u_i||^2 + (ay,y(T))
% s.t.
%       y_t(x,t) - /\y(x,t) + c y^3(x,t) + f(x,t) = sum(i=1:4)au_i(x) u_i(t) 
%       y(x,0) = y_0(x) 
%       dy/dnu(x,t) = 0
%
%       -1 <= u_i(t) <= 1     , i = 1:4
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Numerical methods:    BFGS-method, Finite Element Method,
%                       
%              
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Input parameters:

% parameters

l_max = 10;

k_max = 10;
si = s0;
sigma = 1e-4;
rho = sigma + 0.2;
gamma = 2;
[n,numcontrols] = size(ul);

% function value old
% J_ul = objvalue_pod(ul, ylpod, yqpod, ay, lambda, M, Mt, POD, Mpod);
% 
% % gradient old
% 
% pl = reshape(plpod,r*n,1);
% gl = zeros(n,numcontrols);
% gl(:,1) = lambda*(ul(:,1) + V1pod*pl);
% gl(:,2) = lambda*(ul(:,2) + V2pod*pl);
% gl(:,3) = lambda*(ul(:,3) + V3pod*pl);
% gl(:,4) = lambda*(ul(:,4) + V4pod*pl);
% gl=reshape(gl,n*numcontrols,1);

uk = ul + si*reshape(dk,n,numcontrols);
%Projection on the admissable set
    for j=1:numcontrols
        for i=1:n
            if uk(i,j)<= ua(i,j)
                uk(i,j)=ua(i,j);
            else if uk(i,j)>=ub(i,j)
                    uk(i,j)=ub(i,j);
                end
            end
        end
    end  
% state and adjoint state via POD-DEIM
yk = stateeqn_pod_deim(deltat, y0pod, c, fxyt, au, uk, Mpod , Kpod , POD, U,P, M );
pk = adjointeqn_pod_deim(n,deltat, yk, yqpod, c, ay, Mpod, Kpod, POD, M,U,P);

% value J(y,u)
J_uk = objvalue_pod(uk, yk, yqpod, ay, lambda, M, Mt, POD, Mpod);

% new gradient
pk = reshape(pk,r*n,1);
gk = zeros(n,numcontrols);
gk(:,1) = lambda*(uk(:,1) + V1pod*pk);
gk(:,2) = lambda*(uk(:,2) + V2pod*pk);
gk(:,3) = lambda*(uk(:,3) + V3pod*pk);
gk(:,4) = lambda*(uk(:,4) + V4pod*pk);
gk=reshape(gk,n*numcontrols,1);

% Compute psi(s) and psi'(s) (s. alg  )
psi_si = J_uk-J_ul-sigma*si*gl'*dk;
psi_prime_si = gk'*dk-sigma*gl'*dk;

% A=0;
l=0;
% while (psi(si) < 0 && psi'(si) <= 0)
while ((l < l_max) && (psi_si < 0) && (psi_prime_si <= 0))
% Phase A
% -------------------------------------------------------------------------
%disp(['***                     Line search Phase A: ' num2str(l) ' --- J(f+sd) = ' num2str(J_uk) '    ( s_i = ' num2str(si) ' )']);
    
    si = gamma*si;
    
    uk = ul + si*reshape(dk,n,numcontrols);
    %Projection on the admissable set
    for j=1:numcontrols
        for i=1:n
            if uk(i,j)<= ua(i,j)
                uk(i,j)=ua(i,j);
            else if uk(i,j)>=ub(i,j)
                    uk(i,j)=ub(i,j);
                end
            end
        end
    end
    % state and adjoint state via POD-DEIM
    yk = stateeqn_pod_deim(deltat, y0pod, c, fxyt, au, uk, Mpod , Kpod , POD, U,P, M );
    pk = adjointeqn_pod_deim(n,deltat, yk, yqpod, c, ay, Mpod, Kpod, POD, M,U,P);

    % value J(y,u)
    J_uk = objvalue_pod(uk, yk, yqpod, ay, lambda, M, Mt, POD, Mpod);

    % new gradient
    pk=reshape(pk,r*n,1);
    gk = zeros(n,numcontrols);
    gk(:,1) = lambda*(uk(:,1) + V1pod*pk);
    gk(:,2) = lambda*(uk(:,2) + V2pod*pk);
    gk(:,3) = lambda*(uk(:,3) + V3pod*pk);
    gk(:,4) = lambda*(uk(:,4) + V4pod*pk);
    gk=reshape(gk,n*numcontrols,1);

    % Compute psi(s) and psi'(s) 
    psi_si = J_uk-J_ul-sigma*si*gl'*dk;
    psi_prime_si = gk'*dk-sigma*gl'*dk;

    l=l+1;
    A=1;
end
% if psi(si) < 0 && |psi'(si)| <= (rho-sigma)|phi'(0)|
if (psi_si < 0 && abs(psi_prime_si) <= (rho-sigma)*abs(gl'*dk))
    sk=si;uk = ul + si*reshape(dk,n,numcontrols);
    %Projection on the admissable set
    for j=1:numcontrols
        for i=1:n
            if uk(i,j)<= ua(i,j)
                uk(i,j)=ua(i,j);
            else if uk(i,j)>=ub(i,j)
                    uk(i,j)=ub(i,j);
                end
            end
        end
    end
    
    % state and adjoint state via POD-DEIM
    yk = stateeqn_pod_deim(deltat, y0pod, c, fxyt, au, uk, Mpod , Kpod , POD, U,P, M );
    pk = adjointeqn_pod_deim(n,deltat, yk, yqpod, c, ay, Mpod, Kpod, POD, M,U,P);

    % value J(y,u)
    J_uk = objvalue_pod(uk, yk, yqpod, ay, lambda, M, Mt, POD, Mpod);

    % new gradient
    pk=reshape(pk,r*n,1);
    gk = zeros(n,numcontrols);
    gk(:,1) = lambda*(uk(:,1) + V1pod*pk);
    gk(:,2) = lambda*(uk(:,2) + V2pod*pk);
    gk(:,3) = lambda*(uk(:,3) + V3pod*pk);
    gk(:,4) = lambda*(uk(:,4) + V4pod*pk);
    gk=reshape(gk,n*numcontrols,1);
    
    % Stop --> return si, yk, pk, uk 
    ykpod_deim = yk;
    pkpod_deim = pk;
    
    return;
else
    % if psi(si) < 0 && psi'(si) > 0
    if psi_si < 0 && psi_prime_si > 0
        a=si;
        b=0;
    % if psi(si) >= 0    
    else
        a=0;
        b=si;
    end
    si=1/2*(a+b);
    
    % new stepsize--> compute new solution set
    uk = ul + si*reshape(dk,n,numcontrols);
    %Projection on the admissable set
    for j=1:numcontrols
        for i=1:n
            if uk(i,j)<= ua(i,j)
                uk(i,j)=ua(i,j);
            else if uk(i,j)>=ub(i,j)
                    uk(i,j)=ub(i,j);
                end
            end
        end
    end
    % state and adjoint state via POD-DEIM
    yk = stateeqn_pod_deim(deltat, y0pod, c, fxyt, au, uk, Mpod , Kpod , POD, U,P, M );
    pk = adjointeqn_pod_deim(n,deltat, yk, yqpod, c, ay, Mpod, Kpod, POD, M,U,P);

    % value J(y,u)
    J_uk = objvalue_pod(uk, yk, yqpod, ay, lambda, M, Mt, POD, Mpod);

    % new gradient
    pk=reshape(pk,r*n,1);
    gk = zeros(n,numcontrols);
    gk(:,1) = lambda*(uk(:,1) + V1pod*pk);
    gk(:,2) = lambda*(uk(:,2) + V2pod*pk);
    gk(:,3) = lambda*(uk(:,3) + V3pod*pk);
    gk(:,4) = lambda*(uk(:,4) + V4pod*pk);
    gk=reshape(gk,n*numcontrols,1);

    % Compute psi(s) and psi'(s)
    psi_si = J_uk-J_ul-sigma*si*gl'*dk;
    psi_prime_si = gk'*dk-sigma*gl'*dk;
    
    % psi(a)
    uk_a = ul + a*reshape(dk,n,numcontrols);
    yk_a = stateeqn_pod_deim(deltat, y0pod, c, fxyt, au, uk_a, Mpod , Kpod , POD, U,P, M );
    J_uk_a = objvalue_pod(uk_a, yk_a, yqpod, ay, lambda, M, Mt, POD, Mpod);
    psi_a = J_uk_a-J_ul-sigma*a*gl'*dk;
    
    k=0;
    % while (psi(si) >= psi(ai) || |psi'(si)| > (rho-sigma)|phi'(0)|)
    while ((psi_si >= psi_a || abs(psi_prime_si) > (rho-sigma)*abs(gl'*dk)) && k < k_max )
        %disp(['***                     Line search Phase B: ' num2str(k) ' --- J(f+sd) = ' num2str(J_uk) '    ( s_i = ' num2str(si) ' )']);
%         if A ==1;
%             J_u = J_uk;
%         else
%             J_u = 0;
%         end
        
        % psi(si) >= psi(a)
        if psi_si >= psi_a
            b=si;
        else  
            % psi'(si)(si-ai) < 0
            if psi_prime_si*(si-a) < 0
                a=si;
                % psi(a) new
                uk_a = ul + a*reshape(dk,n,numcontrols);
                yk_a = stateeqn_pod_deim(deltat, y0pod, c, fxyt, au, uk_a, Mpod , Kpod , POD, U,P, M );
                J_uk_a = objvalue_pod(uk_a, yk_a, yqpod, ay, lambda, M, Mt, POD, Mpod);
                psi_a = J_uk_a-J_ul-sigma*a*gl'*dk;
            else
                b=a;
                a=si;
                % psi(a) new
                uk_a = ul + a*reshape(dk,n,numcontrols);
                yk_a = stateeqn_pod_deim(deltat, y0pod, c, fxyt, au, uk_a, Mpod , Kpod , POD, U,P, M );
                J_uk_a = objvalue_pod(uk_a, yk_a, yqpod, ay, lambda, M, Mt, POD, Mpod);
                psi_a = J_uk_a-J_ul-sigma*a*gl'*dk;
            end
        end
        si=1/2*(a+b);
        % new stepsize--> compute new solution set
        uk = ul + si*reshape(dk,n,numcontrols);
        %Projection on the admissable set
        for j=1:numcontrols
            for i=1:n
                if uk(i,j)<= ua(i,j)
                    uk(i,j)=ua(i,j);
                else if uk(i,j)>=ub(i,j)
                        uk(i,j)=ub(i,j);
                    end
                end
            end
        end
        % state and adjoint state via POD-DEIM
        yk = stateeqn_pod_deim(deltat, y0pod, c, fxyt, au, uk, Mpod , Kpod , POD, U,P, M);
        pk = adjointeqn_pod_deim(n,deltat, yk, yqpod, c, ay, Mpod, Kpod, POD, M,U,P);

        % value J(y,u)
        J_uk = objvalue_pod(uk, yk, yqpod, ay, lambda, M, Mt, POD, Mpod);

        % new gradient
        pk=reshape(pk,r*n,1);
        gk = zeros(n,numcontrols);
        gk(:,1) = lambda*(uk(:,1) + V1pod*pk);
        gk(:,2) = lambda*(uk(:,2) + V2pod*pk);
        gk(:,3) = lambda*(uk(:,3) + V3pod*pk);
        gk(:,4) = lambda*(uk(:,4) + V4pod*pk);
        gk=reshape(gk,n*numcontrols,1);

        % Compute psi(s) and psi'(s)
        psi_si = J_uk-J_ul-sigma*si*gl'*dk;
        psi_prime_si = gk'*dk-sigma*gl'*dk;
        
%         if abs(J_u-J_uk) < 1*1e-6
%             disp('break');
%             break
%         end
        k=k+1;
        A=1;
    end
    sk=si;
    
end
% Stop --> return si, yk, pk, uk 
ykpod_deim = yk;
pkpod_deim = reshape(pk,r,n);

    