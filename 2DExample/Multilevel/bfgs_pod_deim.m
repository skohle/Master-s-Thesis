function [yopt_pod,uiopt_pod,popt_pod] = bfgs_pod_deim(n,r,r_deim,m, x1x2,tri,deltat,lambda,c,fxyt,ay,au,ua,ub,ustart,Mpod,Kpod,Mt,M,y0pod,yqpod,POD,U,P,pod_mode,TOL1,TOL2,s0,linesearch_mode,bfgs_mode)
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
% Numerical methods:    BFGS-method, Finite Element Method, POD-DEIM
%                       
%              
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Input parameters:
%
% m             % number of finite elements in [0,L]
% n             % number of time instances in [0,T]
% x1x2, tri     % points, triangles of FEM dicretization
% lambda        % Tikhonov parameter (small)
% yq            % term yq of objective functional
% ay            % term ay of objective functional
% au            % term au of objective functional
% dxyt          % term dxyt on left-hand side of PDE 
% y0            % initial value
% ua            % lower control bound
% ub            % upper control bound
% c             % parameter c in PDE
% ystart        % start value y^0 in SQP method
% ustart        % start value u^0 in SQP method
% pstart        % start value p^0 in SQP method
% M,K           % FE matrices
% V1,V2,V3,V4,v % matrices for projection term
% deltat        % step size in time
% numcontrol   % number of control functions
% Mt  
% cg_mode
% linesearch_mode
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Output:   optimal solution set (uopt, yopt,popt)
% =========================================================================

% Parameters
countermax = 100;

bfgs_mode=2;

[n,numcontrols] = size(ustart);
uk = ustart;
Bk = speye(n*numcontrols,n*numcontrols);

counter = 1;

diff=Inf;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection term from variational inequality
% ui(t) = -1/lambda int_Omega( au_i(x)p(x,t)dx )
% ==> ui = - Vi*p - v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V1
[iau,jau,sau] = find((Mpod*POD'*M*au(:,1))'./lambda);
iV = [];
jV = [];
sV = [];
for k = 1:n
    iV = [iV,k*iau];
    jV = [jV,(k-1)*r+jau];
    sV = [sV,sau];
end
V1pod = sparse(iV,jV,sV,n,n*r);
% V2
[iau,jau,sau] = find((Mpod*POD'*M*au(:,2))'./lambda);
iV = [];
jV = [];
sV = [];
for k = 1:n
    iV = [iV,k*iau];
    jV = [jV,(k-1)*r+jau];
    sV = [sV,sau];
end
V2pod = sparse(iV,jV,sV,n,n*r);
% V3
[iau,jau,sau] = find((Mpod*POD'*M*au(:,3))'./lambda);
iV = [];
jV = [];
sV = [];
for k = 1:n
    iV = [iV,k*iau];
    jV = [jV,(k-1)*r+jau];
    sV = [sV,sau];
end
V3pod = sparse(iV,jV,sV,n,n*r);
% V4
[iau,jau,sau] = find((Mpod*POD'*M*au(:,4))'./lambda);
iV = [];
jV = [];
sV = [];
for k = 1:n
    iV = [iV,k*iau];
    jV = [jV,(k-1)*r+jau];
    sV = [sV,sau];
end

V4pod = sparse(iV,jV,sV,n,n*r);

vpod = sparse(n,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BFGS-Algorithm
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 Step
% Compute one gradient step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Compute state for startcontrol ustart via POD-DEIM
ykpod = stateeqn_pod_deim(deltat, y0pod, c, fxyt, au, uk, Mpod , Kpod , POD, U,P, M );
% Compute the adjoint state via POD_DEIM
pkpod = adjointeqn_pod_deim(n,deltat, ykpod, yqpod, c, ay, Mpod, Kpod, POD, M,U,P);

% pk as vector (POD)
pk=reshape(pkpod,r*n,1);

gk = zeros(n,numcontrols);
gk(:,1) = lambda*(uk(:,1) + V1pod*pk);
gk(:,2) = lambda*(uk(:,2) + V2pod*pk);
gk(:,3) = lambda*(uk(:,3) + V3pod*pk);
gk(:,4) = lambda*(uk(:,4) + V4pod*pk);

gk=reshape(gk,n*numcontrols,1);

% value J(y,u)
J_uk=objvalue_pod(uk, ykpod, yqpod, ay, lambda, M, Mt, POD, Mpod);

    ul=uk;
    gl=gk;
    dk=-gk;
    
    if linesearch_mode==1 || linesearch_mode==2
         %disp('Bisection resp. Armijo is not suitable for BFGS')
         error('myApp:argChk', 'Bisection resp. Armijo is not suitable for BFGS')

     end
    
    if linesearch_mode==3 %Wolfe-Powell
        %disp('***          ---> Wolfe-Powell');
 
           [si,uk,ykpod,pkpod,J_uk,gk]= linesearch_nWP_pod_deim(r,uk,gk,J_uk,s0,dk,y0pod,c,fxyt,au,M,deltat,yqpod,ay,lambda,Mt,ua,ub,Mpod,Kpod,POD,V1pod,V2pod,V3pod,V4pod,U,P, x1x2,r_deim);
    end
    if linesearch_mode==4 %strict Wolfe Powell
        %disp('***          ---> strict Wolfe-Powell');
           [si,uk,ykpod,pkpod,J_uk,gk]= linesearch_WP_pod_deim(r,uk,gk,J_uk,s0,dk,y0pod,c,fxyt,au,M,deltat,yqpod,ay,lambda,Mt,ua,ub,Mpod,Kpod,POD,V1pod,V2pod,V3pod,V4pod,U,P, x1x2,r_deim);
    end
    
    uk=reshape(uk,n*numcontrols,1); 
    ul=reshape(ul,n*numcontrols,1); 
    
    sk=uk-ul;
    fk=gk-gl;
    
    uk=reshape(uk,n,numcontrols);
    
    % Update BFGS matrix
    [iB,jB,sB]=find(Bk-(Bk*(sk*sk')*Bk)/(sk'*Bk*sk)+(fk*fk')/(fk'*sk));
    Bk=sparse(iB,jB,sB);
    
    
%Start loop
while counter< countermax
    % --------------------------------------------------------------
    %disp('***');
    %disp(['***          BFGS-step ' num2str(counter) ' :']);
    % --------------------------------------------------------------
    J_ul = J_uk; % remind J(y_{k-1},u_{k-1})
    ul=uk;
    yl=ykpod;
    gl=gk; % old gradient
    
    
    
    
    normgradientsqu = gk'*gk;
    %disp(['***          ---> Gradient ||gk||^2 = ' num2str(normgradientsqu)]);
    
    %stopping criteria
     %if  abs(J_ul - J_uk) <= 1*1e-6
%     if (normgradientsqu <= TOL2 && abs(J_ul - J_uk) <= TOL1)
%         break;
%      end  
    
   
    % new direction
    %dk=-Bk\gk ;  
    dk=pcg(Bk,-gk);
 
    dl=dk;
    % new stepsize via Wolfe Powell
    % --> new state, control, adjoint state, obj func value and gradient
    
     if linesearch_mode==1 || linesearch_mode==2
         %disp('Bisection resp. Armijo is not suitable for BFGS')
         error('myApp:argChk', 'Bisection resp. Armijo is not suitable for BFGS')

     end
       
     if linesearch_mode==3 %Wolfe-Powell
        %disp('***          ---> Wolfe-Powell');
           [si,uk,ykpod,pkpod,J_uk,gk] = linesearch_nWP_pod_deim(r,uk,gk,J_uk,s0,dk,y0pod,c,fxyt,au,M,deltat,yqpod,ay,lambda,Mt,ua,ub,Mpod,Kpod,POD,V1pod,V2pod,V3pod,V4pod,U,P , x1x2, r_deim);
     end
    if linesearch_mode==4 %strict Wolfe Powell
        %disp('***          ---> strict Wolfe-Powell');
        [si,uk,ykpod,pkpod,J_uk,gk] = linesearch_WP_pod_deim(r,uk,gk,J_uk,s0,dk,y0pod,c,fxyt,au,M,deltat,yqpod,ay,lambda,Mt,ua,ub,Mpod,Kpod,POD,V1pod,V2pod,V3pod,V4pod,U,P , x1x2, r_deim);
    end
    
    
    % -------------- Stopping criteria --------------------

    diffy = error_numerical_Q(yl,ykpod,Mpod,Mt);
    for i = 1:4
        diffu(i) = (ul(:,i)-uk(:,i))'*Mt*(ul(:,i)-uk(:,i));
    end
    diff = diffy + sum(diffu);
    %disp(['***          --->  ||ul-uk||^2+||yl-yk||^2 = ' num2str(diff)]);
    
    if diff < TOL1
        break;
    end
    
    ul=reshape(ul,n*numcontrols,1);
    uk=reshape(uk,n*numcontrols,1); 
    
    sk=uk-ul;
    fk=gk-gl;
    
    if (mod(counter,3)==0)
    
    %if (skfk <= 0)    
        %s0=s0/2;       
    Bl=Bk;    
        switch bfgs_mode
            case 1% Identity
                %disp('***          ---> restart Bk=I');
                Bk=speye(n*numcontrols,n*numcontrols);
            case 2% keep the diagonal elements
                %disp('***          ---> restart Bk=(\)');
                Bk=sparse(n*numcontrols,n*numcontrols);
                for i=1:n*numcontrols
                    Bk(i,i)=Bl(i,i);
                end
            case 3 %keep the tridiagonal elements
                %disp('***          ---> restart Bk=tridiag(Bk)');    
                Bk=sparse(n*numcontrols,n*numcontrols);
                for i=1:n*numcontrols
                    Bk(i,i)=Bl(i,i);
                    if i~=n*numcontrols
                        Bk(i+1,i)=Bl(i+1,i);
                        Bk(i,i+1)=Bl(i,i+1);
                    end                
                end
            case 4
                %disp('***          ---> restart Bk=sum(BK(:,i)');                
                Bk=sparse(n*numcontrols,n*numcontrols);
                for i=1:n*numcontrols
                    Bk(i,i)=sum(abs(Bl(:,i)));
                end
            case 5
                %disp('***          ---> restart Bk=sum(BK(i,:)');                
                Bk=sparse(n*numcontrols,n*numcontrols);
                for i=1:n*numcontrols
                    Bk(i,i)=sum(abs(Bl(i,:)));
                end
            case 6
                %disp('***          ---> restart Bk=max(BK(:,i)');                
                Bk=sparse(n*numcontrols,n*numcontrols);
                for i=1:n*numcontrols
                    Bk(i,i)=max(Bl(:,i));
                end 
        end
    else
        

        % Update BFGS matrix
        Bk=Bk-(Bk*(sk*sk')*Bk)/(sk'*Bk*sk)+(fk*fk')/(fk'*sk);
        
        Bk=0.5.*(Bk+Bk');
    end
 
    uk=reshape(uk,n,numcontrols);
    
    % display objective function values
    %disp(['***          ---> J(yk,uk) = ' num2str(J_uk) '    ( s_i = ' num2str(si) ' )']);
    
    
    counter=counter+1;
    
end
disp(' ');
disp(['Number of iterations: ' num2str(counter)]);
disp(['Stopping criteria  ||ul-uk||^2+||yl-yk||^2 = ' num2str(diff)]);
disp(['Norm of the gradient gk^T*gk = '  num2str(normgradientsqu)]);
disp(['Objective function value J(yk,uk) = ' num2str(J_uk)]);
yopt_pod = POD*ykpod;
uiopt_pod = reshape(uk,n,numcontrols);
popt_pod = POD*reshape(pkpod,r,n);

