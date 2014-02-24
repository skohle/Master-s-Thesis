function [yopt_fem,uiopt_fem,popt_fem] = bfgs_fem(m,deltat,lambda,c,fxyt,ay,au,ua,ub,ustart,M,K,y0,yq,V1,V2,V3,V4,v,Mt,TOL1,TOL2,s0,linesearch_mode,bfgs_mode)
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
%
% m             % number of finite elements in [0,L]
% n             % number of time instances in [0,T]
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
% TOL1          % tolerance obj funct
% TOL2          % tolerance gradient
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Output:   optimal solution set (uopt_fem, yopt_fem,popt_fem)
% =========================================================================

% Parameters

countermax = 100;
[n,numcontrols] = size(ustart);
uk = ustart;
% startmatrix for BFGS matrix
Bk=speye(n*numcontrols,n*numcontrols);
bfgs_mode=2;
counter = 1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BFGS-Algorithm
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %disp('***');
    %disp(['***          BFGS-step ' num2str(counter) ' :']);
    
   
% Compute state for startcontrol ustart
yk = state_equation(n,m,uk,y0,c,fxyt,au,M,K,deltat);
% Compute the adjoint state
pk = adjoint_equation(n,m,yk,deltat,yq,c,ay,M,K);
%Compute gradient 
% pk as vector
pk=reshape(pk,m*n,1);
gk = zeros(n,numcontrols);
gk(:,1) = lambda*(uk(:,1) + V1*pk);
gk(:,2) = lambda*(uk(:,2) + V2*pk);
gk(:,3) = lambda*(uk(:,3) + V3*pk);
gk(:,4) = lambda*(uk(:,4) + V4*pk);
gk=reshape(gk,n*numcontrols,1);
normgradientsqu = gk'*gk;
%disp(['***          ---> Gradient ||gk||^2 = ' num2str(normgradientsqu)]);
 
% value J(y,u)
J_ul=objvalue(uk,yk,yq,ay,lambda,M,Mt);  
ul=uk;
gl=gk;
dk=-gk;

[si,uk,yk,pk,J_uk,gk] = linesearch_nWP(m,uk,gk,J_ul,s0,dk,y0,c,fxyt,au,M,K,deltat,yq,ay,lambda,Mt,ua,ub,V1,V2,V3,V4);   

    uk=reshape(uk,n*numcontrols,1); 
    ul=reshape(ul,n*numcontrols,1); 
    sk=uk-ul;
    fk=gk-gl;
    
    uk=reshape(uk,n,numcontrols);
    
    % Update BFGS matrix
    [iB,jB,sB]=find(Bk-(Bk*(sk*sk')*Bk)/(sk'*Bk*sk)+(fk*fk')/(fk'*sk));
    Bk=sparse(iB,jB,sB);
counter=2;
while counter< countermax
  
    %---------------------------------------------------------------------
    %disp('***');
    %disp(['***          BFGS-step ' num2str(counter) ' :']);
    
    
    
    normgradientsqu = gk'*gk;
    %disp(['***          ---> Gradient ||gk||^2 = ' num2str(normgradientsqu)]);
 
    ul=uk; 
    yl=yk;
    
    dl = dk;
    gl = gk; %old gradient
    J_ul = J_uk; % remind J(y_{k-1},u_{k-1})
    
    
    % new direction
    dk = pcg(Bk,-gk) ;  
    dk=reshape(dk,n,numcontrols);
    
    dk=reshape(dk,n*numcontrols,1);
    
    if linesearch_mode==1 || linesearch_mode==2
         %disp('Bisection resp. Armijo is not suitable for BFGS')
         error('myApp:argChk', 'Bisection resp. Armijo is not suitable for BFGS')

     end
       
    % new stepsize via Wolfe Powell
    % --> new state, control, adjoint state, obj func value and gradient
    if linesearch_mode==3 %Wolfe-Powell
        %disp('***          ---> Wolfe-Powell');
 
           [si,uk,yk,pk,J_uk,gk] = linesearch_nWP(m,uk,gk,J_ul,s0,dk,y0,c,fxyt,au,M,K,deltat,yq,ay,lambda,Mt,ua,ub,V1,V2,V3,V4);   
    end
    if linesearch_mode==4 %strict Wolfe Powell
        %disp('***          ---> srtict Wolfe-Powell');
           [si,uk,yk,pk,J_uk,gk] = linesearch_WP(m,uk,gk,J_ul,s0,dk,y0,c,fxyt,au,M,K,deltat,yq,ay,lambda,Mt,ua,ub,V1,V2,V3,V4);   
    end
    % stopping criteria
    %      if (normgradientsqu <= TOL2 && abs(J_ul - J_uk) <= TOL1) 
%     % if (abs(J_ul - J_uk) <= TOL)
%         break;
%      end  

    diffy = error_numerical_Q(yl,yk,M,Mt);
    for i = 1:4
        diffu(i) = (ul(:,i)-uk(:,i))'*Mt*(ul(:,i)-uk(:,i));
    end
    diff = diffy + sum(diffu);
    %disp(['***          --->  ||ul-uk||^2+||yl-yk||^2 = ' num2str(diff)]);
    
    if diff < TOL1
        break;
    end
    uk=reshape(uk,n*numcontrols,1); 
     
    
         
    % Restart the BFGS-matrix Bk
    if (mod(counter,5)==0)
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
        ul=reshape(ul,n*numcontrols,1); 
        sk=uk-ul;
        fk=gk-gl;

        % Update BFGS matrix
        Bk=Bk-(Bk*(sk*sk')*Bk)/(sk'*Bk*sk)+(fk*fk')/(fk'*sk);
        Bk=0.5.*(Bk+Bk');
    end
    uk=reshape(uk,n,numcontrols);
    
    % display objective function values
    %disp(['***          ---> J(yk,uk) = ' num2str(J_uk) '    ( s_i = ' num2str(si) ' )']);
    
    %    
    counter=counter+1;
%     end
    
    
end
disp(' ');
disp(['Number of iterations: ' num2str(counter)]);
disp(['Stopping criteria  ||ul-uk||^2+||yl-yk||^2 = ' num2str(diff)]);
disp(['Norm of the gradient gk^T*gk = '  num2str(normgradientsqu)]);
disp(['Objective function value J(yk,uk) = ' num2str(J_uk)]);
yopt_fem = yk;
uiopt_fem = reshape(uk,n,numcontrols);
popt_fem = reshape(pk,m,n);

