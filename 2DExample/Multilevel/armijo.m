function [sk,uk,yk,pk,J_uk]= armijo(yk,uk,J_uk,dk,gl,y0,c,fxyt,au,M,K,deltat,yq,ay,lambda,Mt,ua,ub)

% =========================================================================
% Author: Sophia Kohle, Technische Universität Berlin
% =========================================================================
%
[m,n]=size(yk);
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
    
    yk = state_equation(n,m,uk,y0,c,fxyt,au,M,K,deltat);
    
    J_uk = objvalue(uk,yk,yq,ay,lambda,M,Mt);
    
    disp(['***                     Line search: ' num2str(counter) ' --- J(u+sd) = ' num2str(J_uk) '    ( s_i = ' num2str(si) ' )']);    
    
    if J_uk < (J_ul +delta*si*gl'*dk)
        break
    else
        si=(beta1*si+beta2*si)/2;
    end
    counter=counter+1;
end

sk=si;
pk=adjoint_equation(n,m,yk,deltat,yq,c,ay,M,K);
        