
Start the optimization with start.

The following parameters can be chosen.



mmin = 2000;            % minimal number of finite elements in space discretization
n = 500;                % number of time instances
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