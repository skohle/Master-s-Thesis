function [y,ui,p,m,x1,x2,au,M,V1,V2,V3,V4,v,yq,Mt,ay,y0,fxyt,K,tri] =multilevel_bfgs_fem(mmin, n,  mode, linesearch_mode,ustart_ml , t, lambda, c, deltat,T,TOL1,TOL2,s0,bfgs_mode)

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
% -------------------------------------------------------------------------
% 
% Used algorithms:  CG method, finite
%                   element discretization in space, semi-implicit Euler
%                   method in time, proper orthogonal decomposition
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Input parameters:
%              
% mmin          % minimal number of finite elements in space discretization
% n             % number of time instances in [0,T]
% r             % rank of POD basis
% mode          % 0 = Computation by FE method and by POD method
%               % 1 = Computation by FE method
%               % 2 = Computation by POD method
% errormode     % 0 = no computation of analytical a-posteriori-estimate
%               % 1 = computation of analytical a-posteriori-estimate
%
% cg_mode                 % 1: Update with Hestenes-Stiefel
%                         % 2: Update with Flecher-Reeves
%                         % 3: Update with Polak-Riebere
%                         % 4: Update with Hager Zhang
% 
% linesearch_mode         % 1: new search direction via bisection
%                         % 2: new search direction via ...
%                         
% pod_mode                % 1: Computations by POD
%                         % 2: Computations by POD-DEIM
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% OUTPUT:
%

% uifem             - optimal FEM controls u_1, ..., u_4
% yfem              - optimal FEM state
% pfem              - corresponding FEM adjoint state
% uipod             - optimal POD controls u_1, ..., u_4
% ypod              - optimal POD state
% ppod              - corresponding POd adjoint state
% uiopt             - optimal controls u_1, ..., u_4 (constructed)
% yopt              - optimal state (constructed)
% popt              - corresponding adjoint state (constructed)
%
% =========================================================================

generate_data(mmin,n)

load data.mat 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Computation by FE method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(mode == 1 || mode==0)
    
    % optimization by BFGS
    disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
    disp('    Optimization by using FEM');
    disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
    timestart = cputime;

    [yopt_fem,uiopt_fem,popt_fem] = bfgs_fem(m,deltat,lambda,c,fxyt,ay,au,ua,ub,ustart_ml,M,K,y0,yq,V1,V2,V3,V4,v,Mt,TOL1,TOL2,s0,linesearch_mode,bfgs_mode);

    time_fem = cputime - timestart;
    disp(' ');
    disp(['CPU time - FEM optimization: ' num2str(time_fem) 's']);
    disp(' ');

    % sum(i=1:4)au_i(x)ufem_i(t)
    uopt_fem = zeros(m,n);
    for i = 1:n
        for j = 1:m
            for k = 1:numcontrols
                uopt_fem(j,i) = uopt_fem(j,i) + au(j,k)*uiopt_fem(i,k);
            end
        end
    end
    y = yopt_fem;
    ui = uiopt_fem;
    p = popt_fem;

end

if( mode == 1 || mode==0)
    minvalue_fem = objvalue(uiopt_fem,yopt_fem,yq,ay,lambda,M,Mt);
    disp(' ');
    disp(['Objective function value - FEM: ' num2str(minvalue_fem)]);
end



















