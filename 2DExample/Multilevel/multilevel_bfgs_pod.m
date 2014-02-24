function [y,ui,p,m,x1,x2,au,M,V1,V2,V3,V4,v,yq,Mt,ay,y0,fxyt,K,tri] =multilevel_bfgs_pod(mmin, n, r, mode,linesearch_mode,pod_mode, ustart_ml , t, lambda, c, deltat,T,TOL1,TOL2,s0,usnap,it,bfgs_mode)

% =========================================================================
% Author: Sophia Kohle, Technische Universit�t Berlin
% =========================================================================
%
% Masterthesis
% 2D Example_CG[y,ui,p]
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
% -----------------------% upper control bound--------------------------------------------------
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
% 2) Computation by POD method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( mode == 2 || mode == 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2a) Compuation by POD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if pod_mode==0 || pod_mode==1
        
        disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
        disp(['    Optimization by using POD with basis rank r = ' num2str(r)]);
        disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
        
        timestart = cputime;

        % POD basis for snapshots ysnap = G(usnap)
        % --------------------------------------------------------------------
        timestart_podbasis = cputime;
        mode_d=0; 
        [~,POD] = PODbasis(n,m,r,deltat,c,fxyt,au,usnap,M,K,Mt,y0,0,mode_d);
        time_podbasis = cputime - timestart_podbasis;
        
        % Discretization of Omega by POD Galerkin method
        [Mpod,Kpod,y0pod,yqpod] = assem_POD(M,K,y0,yq,POD);
       
        % optimization by CG
        timestart_podopt = cputime;
        
        [yopt_pod,uiopt_pod,popt_pod] = bfgs_pod(r, x1x2,tri,deltat,lambda,c,fxyt,ay,au,ua,ub,ustart_ml,Mpod,Kpod,Mt,M,y0pod,yqpod,POD,linesearch_mode,TOL1,TOL2,s0,bfgs_mode);
        time_podopt = cputime - timestart_podopt;
        time_pod = cputime - timestart;
        disp(' ');
        disp(['CPU time - POD optimization (incl. POD basis): ' num2str(time_pod) 's']);
        disp(['   --> CPU time - POD basis: ' num2str(time_podbasis) 's']);
        disp(['   --> CPU time - POD Optimization: ' num2str(time_podopt) 's']);
    
        
    end
        y = yopt_pod;
        ui = uiopt_pod;
        p = popt_pod;

end


















