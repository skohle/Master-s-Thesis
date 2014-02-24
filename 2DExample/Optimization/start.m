function start

% =========================================================================
% Author: Sophia Kohle, Technische Universit�t Berlin
% =========================================================================
%
% Masterthesis
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Task: Computation of an optimal control for the problem
%
% min 1/2 ||y-yq||^2 + lambda/2 sum(i=1:4) ||u_i||^2 + (ay,y(T))
% s.t.
%       y_t(x,t) - /\y(x,t) + c * y^3(x,t) + f(x,t) = sum(i=1:4)au_i(x) u_i(t) 
%       y(x,0) = y_0(x) 
%       dy/dnu(x,t) = 0
%
%       ua <= u_i(t) <= ub     , i = 1:4
%
% -------------------------------------------------------------------------
%

clear all

mmin = 2000;               % minimal number of finite elements in space discretization
n = 500;                 % number of time instances
r = 3;                  % rank of POD basis 
r_deim= 10;             % number of interpolation points for DEIM
opt_mode=1;             % 0: Optimization with steepest descent
                        % 1: Optimization with non-linear CG
                        % 2: Optimization with BFGS
                        % 3: Optimization with semi-smooth Newton
                        
mode = 2;               % 0: Computation by FE method and by POD method
                        % 1: Computation by FE method
                        % 2: Computation by POD method

pod_mode=0;             % 0: Computation by POD and POD-DEIM
                        % 1: Computations by POD
                        % 2: Computations by POD-DEIM
                                                
errormode = 0;          % 0: no computation of analytical a-posteriori-estimate
                        % 1: computation of analytical a-posteriori-estimate 
                        %     (very complex)
    
linesearch_mode=1;      % 1: new search direction via bisection
                        % 2: new search direction via armijo
                        % 3: new search direction via Wolfe-Powell
                        % 4: new search direction via strict Wolfe-Powell                        
                        
cg_mode=4;              % 1: Update with Hestenes-Stiefel
                        % 2: Update with Flecher-Reeves
                        % 3: Update with Polak-Riebere
                        % 4: Update with Hager Zhang  
                        


switch opt_mode
    case 0
        disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
        disp('  Start optimization with the steepest descent method');
        disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
        example2D_LS(mmin, n, r,r_deim, mode, errormode,linesearch_mode,pod_mode);
    case 1
        disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
        disp('  Start optimization with the non-linear CG method');
        disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
        example2D_CG(mmin, n, r,r_deim, mode, errormode,cg_mode,linesearch_mode,pod_mode);
    case 2
        disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
        disp('  Start optimization with the BFGS method');   
        disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
        example2D_BFGS(mmin, n, r, r_deim, mode, errormode, linesearch_mode, pod_mode);
    case 3
        disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
        disp('  Start optimization with the semi-smooth Newton method');   
        disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
        example2D_ssn(mmin, n, r,r_deim, mode, errormode,pod_mode);
end
   

                  
                        
                        
                        
                        
                        