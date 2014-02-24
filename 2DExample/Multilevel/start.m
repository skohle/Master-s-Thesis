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

mminvec = [180,600,2000];
        %mminvec = [180,600,2000];,10000]   % discretisation in space
nvec = [50,200,500];
        %nvec = [50,200,500,2000];
r = 3;                  % rank of POD basis 
r_deim= 10;             % number of interpolation points for DEIM
mode = 2;               % 0: Computation by FE method and by POD method
                        % 1: Computation by FE method
                        % 2: Computation by POD method
pod_mode=2;             % 0: Computation by POD and POD-DEIM
                        % 1: Computations by POD
                        % 2: Computations by POD-DEIM
opt_mode=0;             % 0: Optimization with steepest descent
                        % 1: Optimization with non-linear CG
                        % 2: Optimization with BFGS                        
                        
errormode = 0;          % 0: no computation of analytical a-posteriori-estimate
                        % 1: computation of analytical a-posteriori-estimate 
                        %     (very complex)
    
linesearch_mode=2;      % 1: new search direction via bisection
                        % 2: new search direction via armijo
                        % 3: new search direction via Wolfe-Powell
                        % 4: new search direction via strict Wolfe-Powell
                                          
                        
cg_mode=4;              % 1: Update with Hestenes-Stiefel
                        % 2: Update with Flecher-Reeves
                        % 3: Update with Polak-Riebere
                        % 4: Update with Hager Zhang  

bfgs_mode=2;
                        
    switch opt_mode
        case 0
            example2D_ml_LS(r,r_deim,mode,errormode,linesearch_mode,pod_mode,mminvec,nvec);
        case 1
            example2D_ml_CG(r,r_deim,mode,errormode,cg_mode,linesearch_mode,pod_mode,mminvec,nvec);
        case 2
            example2D_ml_BFGS(r,r_deim,mode,errormode,cg_mode,linesearch_mode,pod_mode,mminvec,nvec,bfgs_mode);
        case 3

        case 4

    end
           
                        
                        
                        