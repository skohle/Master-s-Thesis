function [sk,uk,ykpod,pkpod,J_uk,gk] = linesearch_nWP_pod(r,ylpod,ul,plpod,J_uk,s0,dk,y0pod,c,fxyt,au,M,deltat,yqpod,ay,lambda,Mt,ua,ub,Mpod,Kpod,POD,Mpod4,V1pod,V2pod,V3pod,V4pod)


% parameters

l_max = 10;

k_max = 5;
si = s0;
sigma = 1e-4;
rho = sigma + 0.2;
gamma = 2;
[n,numcontrols] = size(ul);

% function value old
J_ul = J_uk;

% gradient old

pl = reshape(plpod,r*n,1);
gl = zeros(n,numcontrols);
gl(:,1) = lambda*(ul(:,1) + V1pod*pl);
gl(:,2) = lambda*(ul(:,2) + V2pod*pl);
gl(:,3) = lambda*(ul(:,3) + V3pod*pl);
gl(:,4) = lambda*(ul(:,4) + V4pod*pl);
gl=reshape(gl,n*numcontrols,1);

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
yk = stateeqn_pod(deltat, y0pod, c, fxyt, au, uk, Mpod , Kpod , POD , M , Mpod4);
pk = adjointeqn_pod(n,deltat, yk, yqpod, c, ay, Mpod, Kpod, POD, Mpod4, M);
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

psi_si = J_uk-J_ul-sigma*si*gl'*dk;
%psi_prime_si = gk'*dk-sigma*gl'*dk;
phi_prime_si=gk'*dk;
phi_prime_0=gl'*dk;

l=0;
% while (psi(si) < 0 && psi'(si) <= 0)
while ((l < l_max) && (psi_si <= 0) && (phi_prime_si < rho*phi_prime_0))
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
    yk = stateeqn_pod(deltat, y0pod, c, fxyt, au, uk, Mpod , Kpod , POD , M , Mpod4);
    pk = adjointeqn_pod(n,deltat, yk, yqpod, c, ay, Mpod, Kpod, POD, Mpod4, M);

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

    psi_si = J_uk-J_ul-sigma*si*gl'*dk;
    phi_prime_si=gk'*dk;
    phi_prime_0=gl'*dk;
    l=l+1;
end
% if psi(si) < 0 && |psi'(si)| <= (rho-sigma)|phi'(0)|
if (psi_si <= 0)
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
    yk = stateeqn_pod(deltat, y0pod, c, fxyt, au, uk, Mpod , Kpod , POD , M , Mpod4);
    pk = adjointeqn_pod(n,deltat, yk, yqpod, c, ay, Mpod, Kpod, POD, Mpod4, M);

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
    
    ykpod = yk;
    pkpod = pk;
    %disp(['***  --- J(f+sd) = ' num2str(J_uk) '    ( s_i = ' num2str(si) ' )']);
  
    return;
else
    a=0;
    b=si;

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
    yk = stateeqn_pod(deltat, y0pod, c, fxyt, au, uk, Mpod , Kpod , POD , M , Mpod4);
    pk = adjointeqn_pod(n,deltat, yk, yqpod, c, ay, Mpod, Kpod, POD, Mpod4, M);

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

    psi_si = J_uk-J_ul-sigma*si*gl'*dk;
    phi_prime_si=gk'*dk;
    phi_prime_0=gl'*dk;
    
    k=0;
    % while (psi(si) >= psi(ai) || |psi'(si)| > (rho-sigma)|phi'(0)|)
    while (((psi_si > 0) || (phi_prime_si < (rho*phi_prime_0))) && k < k_max )
      %  disp(['***                     Line search Phase B: ' num2str(k) ' --- J(f+sd) = ' num2str(J_uk) '    ( s_i = ' num2str(si) ' )']);

        % psi(si) >= psi(a)
        if psi_si >0
            b=si;
        else  
            a=si;
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
        yk = stateeqn_pod(deltat, y0pod, c, fxyt, au, uk, Mpod , Kpod , POD , M , Mpod4);
        pk = adjointeqn_pod(n,deltat, yk, yqpod, c, ay, Mpod, Kpod, POD, Mpod4, M);

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

        psi_si = J_uk-J_ul-sigma*si*gl'*dk;
        phi_prime_si=gk'*dk;
        phi_prime_0=gl'*dk;
        
        k=k+1;
    end
    sk=si;
    %disp(['***  --- J(f+sd) = ' num2str(J_uk) '    ( s_i = ' num2str(si) ' )']);
end

ykpod = yk;
pkpod = reshape(pk,r,n);

    