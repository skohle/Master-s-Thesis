function [sk,uk,yk,pk,J_uk]= armijo_pod(yk,uk,J_uk,dk,gl,y0pod,c,fxyt,au,M,deltat,yqpod,ay,lambda,Mt,ua,ub, Mpod , Kpod , POD ,  Mpod4)

% =========================================================================
% Author: Sophia Kohle, Technische Universität Berlin
% =========================================================================
%

% Follow the Armijo Algorithm described in chapter 3.2.2

[n,numcontrol]=size(uk);
delta=0.001;

si=-(gl'*dk)/(dk'*dk);

beta1=0.1;
beta2=0.09;
countermax=5;
counter=1;

while counter < countermax
    J_ul=J_uk;
    ul=uk;
    yl=yk;
    uk = ul+si*reshape(dk,n,numcontrol);
    %Projection on the admissable set
    for j=1:numcontrol
        for k=1:n
            if uk(k,j)<= ua(k,j)
                uk(k,j)=ua(k,j);
            else if uk(k,j)>=ub(k,j)
                    uk(k,j)=ub(k,j);
                end
            end
        end
    end
    
    yk = stateeqn_pod(deltat, y0pod, c, fxyt, au, uk, Mpod , Kpod , POD , M , Mpod4);
    
    J_uk = objvalue_pod(uk, yk, yqpod, ay, lambda, M, Mt, POD, Mpod);
    
    disp(['***                     Line search: ' num2str(counter) ' --- J(u+sd) = ' num2str(J_uk) '    ( s_i = ' num2str(si) ' )']);    
    % stopping criteria
    if J_uk < (J_ul +delta*si*gl'*dk)
        break
    else
        si=(beta1*si+beta2*si)/2;
    end
    counter=counter+1;
end

sk=si;
pk=adjointeqn_pod(n,deltat, yk, yqpod, c, ay, Mpod, Kpod, POD, Mpod4, M);
        