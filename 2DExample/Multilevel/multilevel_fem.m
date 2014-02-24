function [y,ui,p,m,x1,x2,au,M,V1,V2,V3,V4,v,yq,Mt,ay,y0,fxyt,K,tri] =multilevel_bfgs_fem(mmin, n, r,r_deim, mode, errormode,linesearch_mode,pod_mode, ustart , t, lambda, c, deltat,T,TOL1,TOL2,s0)

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

omegashape = rectangle_g(0,pi,0,pi);          % geometry file for space Omega
gammashape = rectangle_b('0', '0', '0', '0'); % boundary file for space Omega

[n,numcontrols]=size(ustart);

ua = -ones(n,numcontrols);   % lower control bound
ub = ones(n,numcontrols);    % upper control bound
usnap(:,1) = 0.5*ones(n,1);  % controls for generation of snapshots
usnap(:,2) = -0.5*ones(n,1);
usnap(:,3) = 0.5*ones(n,1);
usnap(:,4) = -0.5*ones(n,1);


% Discretization of (0,T) piecewise constant

Mt = deltat * eye(n,n);
Mt(1,1) = Mt(1,1)/2;
Mt(n,n) = Mt(n,n)/2;
   
% Discretization of Omega by Finite Element method
[m,x1x2,edges,tri] = grid_generation(mmin, omegashape);
disp(' ');
disp(['   Number of Finite Elements: m = ' num2str(m)]); 
disp(' ');

x1 = x1x2(1,:)';
x2 = x1x2(2,:)';

au(:,1) = max(0,10-50*(x1-pi/4).^2-50*(x2-pi/4).^2);      % term au of objective functional
au(:,2) = max(0,10-50*(x1-3/4*pi).^2-50*(x2-pi/4).^2);
au(:,3) = max(0,10-50*(x1-3/4*pi).^2-50*(x2-3/4*pi).^2);
au(:,4) = max(0,10-50*(x1-pi/4).^2-50*(x2-3/4*pi).^2);
[M,K,~,~,~,~,~] = assem_FEM(x1x2, edges, tri, gammashape, au);

% Parameters (2)
% ------------------------------------------------------------------------
ay = cos(x1).*cos(x2)*lambda*T^2;	% term ay of objective functional
y0 = cos(x1).*cos(x2);            % initial value      

yq = zeros(m,n);                % term yq of objective functional
for i = 1:m
    for j = 1:n
        yq(i,j) = cos(x1(i))*cos(x2(i))*(1 + 2*lambda*t(j) - 2*lambda*t(j)^2 - 3*c*lambda*t(j)^2*cos(x1(i))^2*cos(x2(i))^2);
    end
end

% Parameters (3)
pxy = cos(x1).*cos(x2);
        Aaup = au'*M*pxy;
    rei  = zeros(n,numcontrols);
    for i = 1:numcontrols
        rei(:,i) = min(ub(:,i),max(ua(:,i),-t.^2*Aaup(i)));
    end

        re = zeros(m,n);
        for i = 1:n
            for j = 1:m
                for k = 1:numcontrols
                    re(j,i) = re(j,i) + au(j,k)*rei(i,k);
                end
            end
        end


fxyt = zeros(m,n);              % term f(x,y,t) on left-hand side of PDE 
for i = 1:n
    fxyt(:,i) = re(:,i) - 2*cos(x1).*cos(x2) - c*cos(x1).^3.*cos(x2).^3;
end
fxyt = sparse(fxyt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection term from variational inequality
% ui(t) = -1/lambda int_Omega( au_i(x)p(x,t)dx )
% ==> ui = - Vi*p - v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V1
[iF,jF,sF] = find((au(:,1)'*M)./lambda);
iV = [];
jV = [];
sV = [];
for k = 1:n
    iV = [iV,k*iF];
    jV = [jV,(k-1)*m+jF];
    sV = [sV,sF];
end
V1 = sparse(iV,jV,sV,n,n*m);
% V2
[iF,jF,sF] = find((M*au(:,2))'./lambda);
iV = [];
jV = [];
sV = [];
for k = 1:n
    iV = [iV,k*iF];
    jV = [jV,(k-1)*m+jF];
    sV = [sV,sF];
end
V2 = sparse(iV,jV,sV,n,n*m);
% V3
[iF,jF,sF] = find((M*au(:,3))'./lambda);
iV = [];
jV = [];
sV = [];
for k = 1:n
    iV = [iV,k*iF];
    jV = [jV,(k-1)*m+jF];
    sV = [sV,sF];
end
V3 = sparse(iV,jV,sV,n,n*m);
% V4
[iF,jF,sF] = find((M*au(:,4))'./lambda);
iV = [];
jV = [];
sV = [];
for k = 1:n
    iV = [iV,k*iF];
    jV = [jV,(k-1)*m+jF];
    sV = [sV,sF];
end
V4 = sparse(iV,jV,sV,n,n*m);

v = sparse(n,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Computation by FE method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(mode == 1 || mode==0)
    
    % optimization by BFGS
    disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
    disp('    Optimization by using FEM');
    disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
    timestart = cputime;

    [yopt_fem,uiopt_fem,popt_fem] = bfgs_fem(m,deltat,lambda,c,fxyt,ay,au,ua,ub,ustart,M,K,y0,yq,V1,V2,V3,V4,v,Mt,TOL1,TOL2,s0,linesearch_mode);

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



















