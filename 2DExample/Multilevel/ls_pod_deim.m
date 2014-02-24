function [yopt_pod,uopt_pod,popt_pod] = ls_pod_deim(n,r,r_deim, deltat,lambda,c,fxyt,ay,au,ua,ub,ustart,Mpod,Kpod,Mt,M,y0pod,yqpod,POD,U,P,linesearch_mode,TOL1,TOL2,sk)
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
% Numerical methods:    CG-method, POD-DEIM
%                       
%              
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Input parameters:
%
% m             % number of finite elements in [0,L]
% n             % number of time instances in [0,T]
% x1x2, tri       % points, triangles of FEM dicretization
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
% numcontrols   % number of control functions
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
% Parameters
countermax = 100;

numcontrols = size(au,2);

uk = ustart;

counter = 1;

reset= 5;


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start CG
% 1. Solve the non-linear state equation  
% 2. Solve the non-linear adjoint equation  
% 3. Compute the gradient and the update
% 4. Compute a new stepsize and the next serach direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%State and adjoint state in rxr!!!
ykpod = stateeqn_pod_deim(deltat, y0pod, c, fxyt, au, uk, Mpod , Kpod , POD, U,P, M);

pkpod = adjointeqn_pod_deim(n,deltat, ykpod, yqpod, c, ay, Mpod, Kpod, POD, M,U,P);

J_uk=objvalue_pod(uk, ykpod, yqpod, ay, lambda, M, Mt, POD, Mpod);

while counter < countermax

    J_ul = J_uk; % old

    %disp('***');
    %disp(['***          Steepest descent-step ' num2str(counter) ' :']);
    
        pk=reshape(pkpod,r*n,1);

        gk = zeros(n,numcontrols);
        gk(:,1) = lambda*(uk(:,1) + V1pod*pk);
        gk(:,2) = lambda*(uk(:,2) + V2pod*pk);
        gk(:,3) = lambda*(uk(:,3) + V3pod*pk);
        gk(:,4) = lambda*(uk(:,4) + V4pod*pk);
        gk=reshape(gk,n*numcontrols,1);
        dk=-gk;%/norm(gk);
        
        gl=gk; %Old Gradient
        dl=dk; %Old 
        yl=ykpod;
        ul=uk;
    
   normgradientsqu = gk'*gk;
    %disp(['***          ---> Gradient ||gk||^2 = ' num2str(normgradientsqu)]);
    %initial stepsize
    %sk=sk*(gl'*dl)/(gk'*dk);
    switch linesearch_mode
        case 1
            if mod(counter,reset)==0
                sk=sk/2;
            end
            [~,uk,ykpod,pkpod,J_uk]= linesearch_half_pod_deim(ykpod,uk,J_uk,sk,dk,y0pod,c,fxyt,au,M,deltat,yqpod,ay,lambda,Mt,ua,ub,Mpod,Kpod,POD,P,U);
        case 2
            %disp('***          Armijo stepsize rule');
            if counter == 1
                [~,uk,ykpod,pkpod,J_uk]= linesearch_half_pod_deim(ykpod,uk,J_uk,sk,dk,y0pod,c,fxyt,au,M,deltat,yqpod,ay,lambda,Mt,ua,ub,Mpod,Kpod,POD,P,U);
            else
                [sk,uk,ykpod,pkpod,J_uk]= armijo_pod_deim(ykpod,uk,J_uk,dk,gl,c,fxyt,au,M,deltat,yqpod,ay,lambda,Mt,ua,ub,y0pod, Mpod , Kpod , POD ,  U,P);           
            end
        case 3 %Wolfe-Powell
            %disp('***          ---> Wolfe-Powell');
            [si,uk,ykpod,pkpod,J_uk,gk]= linesearch_nWP_pod_deim(r,uk,gk,J_uk,dk,gl,y0pod,c,fxyt,au,M,deltat,yqpod,ay,lambda,Mt,ua,ub,Mpod,Kpod,POD,V1pod,V2pod,V3pod,V4pod,U,P, x1x2,r_deim);
        case 4 %strict Wolfe Powell
            %disp('***          ---> strict Wolfe-Powell');
            [si,uk,ykpod,pkpod,J_uk,gk]= linesearch_WP_pod_deim(r,uk,gk,J_uk,s0,dk,y0pod,c,fxyt,au,M,deltat,yqpod,ay,lambda,Mt,ua,ub,Mpod,Kpod,POD,V1pod,V2pod,V3pod,V4pod,U,P, x1x2,r_deim);
    end

% display objective function values
%disp(['***          ---> J(yk,uk) = ' num2str(J_uk)]);

%  if (normgradientsqu <= TOL2 && abs(J_ul - J_uk) <= TOL1) 
% % if (abs(J_ul - J_uk) <= TOL)
%     break;
%  end     

% stopping condition
    diffy = error_numerical_Q(yl,ykpod,Mpod,Mt);
    for i = 1:4
        diffu(i) = (ul(:,i)-uk(:,i))'*Mt*(ul(:,i)-uk(:,i));
    end
    diff = diffy + sum(diffu);
    %disp(['***          --->  ||ul-uk||^2+||yl-yk||^2 = ' num2str(diff)]);

    % if ||ul-uk||^2+||yl-yk||^2 <= TOL1
    if diff < TOL1 && normgradientsqu<= TOL2
        break;
    end
counter=counter+1;  
end
disp(' ');
disp(['Number of iterations: ' num2str(counter)]);
disp(['Stopping criteria  ||ul-uk||^2+||yl-yk||^2 = ' num2str(diff)]);
disp(['Norm of the gradient gk^T*gk = '  num2str(normgradientsqu)]);
disp(['Objective function value J(yk,uk) = ' num2str(J_uk)]);



yopt_pod = POD*ykpod;
uopt_pod = uk;
popt_pod = POD*pkpod;

