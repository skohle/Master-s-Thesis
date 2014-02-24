function [yopt,uopt,popt] = ls_fem(m,deltat,lambda,c,fxyt,ay,au,ua,ub,ustart,M,K,y0,yq,V1,V2,V3,V4,v,Mt,linesearch_mode,TOL1,TOL2,sk)

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
% Numerical methods:    Steepest descent-method, Finite Element Method,
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
% linesearch_mode
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Output:   optimal solution set (uopt, yopt,popt)
% =========================================================================

% Parameters

countermax = 100;

uk = ustart;
[n,numcontrols] = size(uk);
diff=Inf;

counter = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start steepest decent
% 1. Solve the non-linear state equation  
% 2. Solve the non-linear adjoint equation  
% 3. Compute the gradient 
% 4. Compute a new stepsize and the next search direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Solve the non-linear state equation  
yk = state_equation(n,m,uk,y0,c,fxyt,au,M,K,deltat);
% 2. Solve the non-linear adjoint equation  
pk = adjoint_equation(n,m,yk,deltat,yq,c,ay,M,K);
% objective value J(y,u)
J_uk=objvalue(uk,yk,yq,ay,lambda,M,Mt);

% Loop
while counter < countermax
    
    %disp('***');
    %disp(['***          Steepest decent-step ' num2str(counter) ' :']);
  
    %Store data of the previous step
    J_ul = J_uk; 
    ul=uk; 
    pl=pk; 
    yl=yk; 
    if counter ~= 1
        gl=gk;
        dl=dk; 
    end
    
    % Compute the gradient
    pk=reshape(pk,m*n,1);
    gk = zeros(n,numcontrols);
    gk(:,1) = lambda*(uk(:,1) + V1*pk);
    gk(:,2) = lambda*(uk(:,2) + V2*pk);
    gk(:,3) = lambda*(uk(:,3) + V3*pk);
    gk(:,4) = lambda*(uk(:,4) + V4*pk);
    dk=-gk;
    gk=reshape(gk,n*numcontrols,1);

    dk=reshape(dk,n*numcontrols,1); 
        
    
    % euclidean norm of the gradient
    normgradientsqu = gk'*gk;
    %disp(['***          ---> Gradient ||gk||^2 = ' num2str(normgradientsqu)]);
    
    % Choose stepsize rule
    switch linesearch_mode      
        case 1 %bisection
            %disp('***               Bisection stepsize rule');
            %initial stepsize
            if counter ~= 1
                sk=sk*(gl'*dl)/(gk'*dk);
            end
            [sk,uk,yk,pk,J_uk]= linesearch_half(yk,uk,J_uk,sk,dk,y0,c,fxyt,au,M,K,deltat,yq,ay,lambda,Mt,ua,ub);
        case 2 %Armijo
            %disp('***                Armijo stepsize rule');
            if counter == 1 
                [sk,uk,yk,pk,J_uk]= linesearch_half(yk,uk,J_uk,sk,dk,y0,c,fxyt,au,M,K,deltat,yq,ay,lambda,Mt,ua,ub);
            else
                [sk,uk,yk,pk,J_uk]= armijo(yk,uk,J_uk,dk,gl,y0,c,fxyt,au,M,K,deltat,yq,ay,lambda,Mt,ua,ub);
            end
        case 3 %Wolfe-Powell
            %disp('***                Wolfe-Powell');
            [si,uk,yk,pk,J_uk,gk] = linesearch_nWP(m,uk,gk,J_ul,sk,dk,y0,c,fxyt,au,M,K,deltat,yq,ay,lambda,Mt,ua,ub,V1,V2,V3,V4);   
        case 4 %strict Wolfe Powell
            %disp('***                strict Wolfe-Powell');
            [si,uk,yk,pk,J_uk,gk] = linesearch_WP(m,uk,gk,J_ul,sk,dk,y0,c,fxyt,au,M,K,deltat,yq,ay,lambda,Mt,ua,ub,V1,V2,V3,V4);   
    end
    % display objective function values
    %disp(['***          ---> J(yk,uk) = ' num2str(J_uk)]);
    
    
    % -------------- Stopping criteria --------------------
%      if (normgradientsqu <= TOL2 && abs(J_ul - J_uk) <= TOL1) 
%    %  if (abs(J_ul - J_uk) <= TOL)
%         break;
%      end    
    
    diffy = error_numerical_Q(yl,yk,M,Mt);
    for i = 1:4
        diffu(i) = (ul(:,i)-uk(:,i))'*Mt*(ul(:,i)-uk(:,i));
    end
    diff = diffy + sum(diffu);
    %disp(['***          --->  ||ul-uk||^2+||yl-yk||^2 = ' num2str(diff)]);
    
    % if ||ul-uk||^2+||yl-yk||^2 <= TOL1
    if diff <= TOL2
        break;
    end    
    counter=counter+1;  
end
disp(' ');
disp(['Number of iterations: ' num2str(counter)]);
disp(['Stopping criteria  ||ul-uk||^2+||yl-yk||^2 = ' num2str(diff)]);
disp(['Norm of the gradient gk^T*gk = '  num2str(normgradientsqu)]);
disp(['Objective function value J(yk,uk) = ' num2str(J_uk)]);
yopt = yk;
uopt = uk;
popt = pk;
