function [sk,uk,yk,pk,J_uk,gk] = linesearch_WP(m,ul,gl,J_ul,s0,dk,y0,c,fxyt,au,M,K,deltat,yq,ay,lambda,Mt,ua,ub,V1,V2,V3,V4)
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

k_max = 5;
si = s0;
sigma = 1e-4;
rho = sigma + 0.01;
gamma = 1.5;
[n,numcontrols] = size(ul);

% new control 
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
% state and adjoint state via FEM
yk = state_equation(n,m,uk,y0,c,fxyt,au,M,K,deltat);
pk = adjoint_equation(n,m,yk,deltat,yq,c,ay,M,K);

% value J(y,u)
J_uk=objvalue(uk,yk,yq,ay,lambda,M,Mt);

% new gradient
pk=reshape(pk,m*n,1);
gk = zeros(n,numcontrols);
gk(:,1) = lambda*(uk(:,1) + V1*pk);
gk(:,2) = lambda*(uk(:,2) + V2*pk);
gk(:,3) = lambda*(uk(:,3) + V3*pk);
gk(:,4) = lambda*(uk(:,4) + V4*pk);
gk=reshape(gk,n*numcontrols,1);

% Compute psi(s) and psi'(s) (s. alg  )
psi_si = J_uk-J_ul-sigma*si*gl'*dk;
psi_prime_si = gk'*dk-sigma*gl'*dk;

A=0;

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
    % state and adjoint state via FEM
    yk = state_equation(n,m,uk,y0,c,fxyt,au,M,K,deltat);
    pk = adjoint_equation(n,m,yk,deltat,yq,c,ay,M,K);
    % value J(y,u)
    J_uk=objvalue(uk,yk,yq,ay,lambda,M,Mt);

    % new gradient
    pk=reshape(pk,m*n,1);
    gk = zeros(n,numcontrols);
    gk(:,1) = lambda*(uk(:,1) + V1*pk);
    gk(:,2) = lambda*(uk(:,2) + V2*pk);
    gk(:,3) = lambda*(uk(:,3) + V3*pk);
    gk(:,4) = lambda*(uk(:,4) + V4*pk);
    gk=reshape(gk,n*numcontrols,1);

    % Compute psi(s) and psi'(s)
    psi_si = J_uk-J_ul-sigma*si*gl'*dk;
    psi_prime_si = gk'*dk-sigma*gl'*dk;
    
    l=l+1;
    A=1;
end
% if psi(si) < 0 && |psi'(si)| <= (rho-sigma)|phi'(0)|
if (psi_si < 0 && abs(psi_prime_si) <= (rho-sigma)*abs(gl'*dk))
    sk=si;
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
    % state and adjoint state via FEM
    yk = state_equation(n,m,uk,y0,c,fxyt,au,M,K,deltat);
    pk = adjoint_equation(n,m,yk,deltat,yq,c,ay,M,K);
    % value J(y,u)
    J_uk=objvalue(uk,yk,yq,ay,lambda,M,Mt);
    
    % new gradient
    pk=reshape(pk,m*n,1);
    gk = zeros(n,numcontrols);
    gk(:,1) = lambda*(uk(:,1) + V1*pk);
    gk(:,2) = lambda*(uk(:,2) + V2*pk);
    gk(:,3) = lambda*(uk(:,3) + V3*pk);
    gk(:,4) = lambda*(uk(:,4) + V4*pk);
    gk=reshape(gk,n*numcontrols,1);
%     
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
    % state and adjoint state via FEM
    yk = state_equation(n,m,uk,y0,c,fxyt,au,M,K,deltat);
    pk = adjoint_equation(n,m,yk,deltat,yq,c,ay,M,K);
    
    % value J(y,u)
    J_uk=objvalue(uk,yk,yq,ay,lambda,M,Mt);

    % new gradient
    pk=reshape(pk,m*n,1);
    gk = zeros(n,numcontrols);
    gk(:,1) = lambda*(uk(:,1) + V1*pk);
    gk(:,2) = lambda*(uk(:,2) + V2*pk);
    gk(:,3) = lambda*(uk(:,3) + V3*pk);
    gk(:,4) = lambda*(uk(:,4) + V4*pk);
    gk=reshape(gk,n*numcontrols,1);

    % Compute psi(s) and psi'(s)
    psi_si = J_uk-J_ul-sigma*si*gl'*dk;
    psi_prime_si = gk'*dk-sigma*gl'*dk;
    
    % psi(a)
    uk_a = ul + a*reshape(dk,n,numcontrols);
    yk_a = state_equation(n,m,uk_a,y0,c,fxyt,au,M,K,deltat);
    J_uk_a=objvalue(uk_a,yk_a,yq,ay,lambda,M,Mt);
    psi_a = J_uk_a-J_ul-sigma*a*gl'*dk;
    
    k=0;
    % while (psi(si) >= psi(ai) || |psi'(si)| > (rho-sigma)|phi'(0)|)
    while ((psi_si >= psi_a || abs(psi_prime_si) > (rho-sigma)*abs(gl'*dk)) && k < k_max )
        if A ==1;
            J_u = J_uk;
        else
            J_u = 0;
        end
        
        %disp(['***                     Line search Phase B: ' num2str(k) ' --- J(f+sd) = ' num2str(J_uk) '    ( s_i = ' num2str(si) ' )']);

        % psi(si) >= psi(a)
        if psi_si >= psi_a
            b=si;
        else  
            % psi'(si)(si-ai) < 0
            if psi_prime_si*(si-a) < 0
                a=si;
                % psi(a) new
                uk_a = ul + a*reshape(dk,n,numcontrols);
                yk_a = state_equation(n,m,uk_a,y0,c,fxyt,au,M,K,deltat);
                J_uk_a=objvalue(uk_a,yk_a,yq,ay,lambda,M,Mt);
                psi_a = J_uk_a-J_ul-sigma*a*gl'*dk;
            else
                b=a;
                a=si;
                % psi(a) new
                uk_a = ul + a*reshape(dk,n,numcontrols);
                yk_a = state_equation(n,m,uk_a,y0,c,fxyt,au,M,K,deltat);
                J_uk_a=objvalue(uk_a,yk_a,yq,ay,lambda,M,Mt);
                psi_a = J_uk_a-J_ul-sigma*a*gl'*dk;
            end
        end
        si=1/2*(a+b);
        % new stepsize--> compute new solution set
        uk = ul + si*reshape(dk,n,numcontrols);
        % Projection on the admissable set
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
        % state and adjoint state via FEM
        yk = state_equation(n,m,uk,y0,c,fxyt,au,M,K,deltat);
        pk = adjoint_equation(n,m,yk,deltat,yq,c,ay,M,K);
        
        % value J(y,u)
        J_uk=objvalue(uk,yk,yq,ay,lambda,M,Mt);

        % new gradient
        pk=reshape(pk,m*n,1);
        gk = zeros(n,numcontrols);
        gk(:,1) = lambda*(uk(:,1) + V1*pk);
        gk(:,2) = lambda*(uk(:,2) + V2*pk);
        gk(:,3) = lambda*(uk(:,3) + V3*pk);
        gk(:,4) = lambda*(uk(:,4) + V4*pk);
        gk=reshape(gk,n*numcontrols,1);
% 
        % Compute psi(s) and psi'(s)
        psi_si = J_uk-J_ul-sigma*si*gl'*dk;
        psi_prime_si = gk'*dk-sigma*gl'*dk;
        
%         if abs(J_ul-J_uk) < 1*1e-6
%             disp('break');
%             break
%         end
        k=k+1;
        A=1;
    end
    sk=si;
    
end

    
    

