function [sk,uk_new,yk_new,pk_new,J_uk_new]= linesearch_half_pod_deim(ykpod,uk,J_uk,s0,dk,y0pod,c,fxyt,au,M,deltat,yqpod,ay,lambda,Mt,ua,ub,Mpod,Kpod,POD,P,U)

% =========================================================================
% Author: Sophia Kohle, Technische Universität Berlin
% =========================================================================
%
% Line search algorithm by using s = s/2, if J(uk,yk) < J(uk+1,yk+1)
% POD
%

%parameters
counter=1;
countermax=10;

[n,numcontrol]=size(uk);

% old step
% yk, uk, pk
    if J_uk == 0 || s0 == 0
        sk = 0;
        uk_new = uk;
        yk_new = ykpod;
        pk_new = pk;
        J_uk_new = J_uk;
        return;
    end
si=s0*2;

J_uk_new = J_uk;

while J_uk <= J_uk_new && J_uk > 0 && counter < countermax   
    si=si/2;
    uk_new = uk+si*reshape(dk,n,numcontrol);
    %Projection on th admissable set
    for j=1:numcontrol
        for i=1:n
            if uk_new(i,j)<= ua(i,j)
                uk_new(i,j)=ua(i,j);
            else if uk_new(i,j)>=ub(i,j)
                    uk_new(i,j)=ub(i,j);
                end
            end
        end
    end 
    
    yk_new = stateeqn_pod_deim(deltat, y0pod, c, fxyt, au, uk_new, Mpod , Kpod , POD, U,P, M );
    J_uk_new = objvalue_pod(uk_new, yk_new, yqpod, ay, lambda, M, Mt, POD, Mpod);
    %disp(['***                     Line search: ' num2str(counter) ' --- J(u+sd) = ' num2str(J_uk_new) '    ( s_i = ' num2str(si) ' )']);    
    counter=counter+1;
end
sk=si;
pk_new=adjointeqn_pod_deim(n,deltat, yk_new, yqpod, c, ay, Mpod, Kpod, POD, M,U,P);