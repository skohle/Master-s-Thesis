function example2D_CG(mmin,n,r,r_deim,mode,errormode,cg_mode,linesearch_mode,pod_mode)

% =========================================================================
% Author: Sophia Kohle, Technische Universit�t Berlin
% =========================================================================
%
% Masterthesis
% 2D Example_CG
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
reset=5;
disp(' ');
switch linesearch_mode
    case 1
        disp('Stepsize search with Bisection stepsize rule');
    case 2
        disp('Stepsize search with Armijo stepsize rule');
    case 3 %Wolfe-Powell
        disp('Stepsize search with Wolfe-Powell');
    case 4 %strict Wolfe Powell
        disp('Stepsize search with strict Wolfe-Powell');
end

switch cg_mode        
    case 1 %Hestenes-Stiefel
        disp('Update formula: Hestenes-Stiefel');
    case 2 % Fletcher-Reeves
        disp('Update formula: Fletcher-Reeves');
    case 3 %Polak-Ribiere
        disp('Update formula: Polak-Ribiere');
    case 4 %Hager-Zhang
        disp('Update formula: Hager-Zhang');
end
disp(' ');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Computation by FE method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(mode == 0 || mode == 1)
    
    % optimization by CG
    disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
    disp('    Optimization by using FEM');
    disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
    timestart = cputime;
    %s0=5;
    TOL2=3.6*1e-3;
    [yopt_fem,uiopt_fem,popt_fem] = cg_fem(m,deltat,lambda,c,fxyt,ay,au,ua,ub,ustart,M,K,y0,yq,V1,V2,V3,V4,v,Mt,cg_mode,linesearch_mode,TOL1,TOL2,s0,reset);

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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Computation by POD method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(mode == 0 || mode == 2)
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
        TOL2=3.14*1e-3;
        [yopt_pod,uiopt_pod,popt_pod] = cg_pod(n,r, x1x2,tri,deltat,lambda,c,fxyt,ay,au,ua,ub,ustart,Mpod,Kpod,Mt,M,y0pod,yqpod,POD,linesearch_mode,cg_mode,TOL1,TOL2,s0,reset);
        
        time_podopt = cputime - timestart_podopt;
        time_pod = cputime - timestart;
        disp(' ');
        disp(['CPU time - POD optimization (incl. POD basis): ' num2str(time_pod) 's']);
        disp(['   --> CPU time - POD basis: ' num2str(time_podbasis) 's']);
        disp(['   --> CPU time - POD Optimization: ' num2str(time_podopt) 's']);
    
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2b) Computation by POD-DEIM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if pod_mode==0 || pod_mode==2
            disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
            disp(['    Optimization by using POD-DEIM with basis rank r = ' num2str(r) ' and POD-DEIM rank=' num2str(r_deim)]); 
            disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
            
            timestart = cputime;
            
            % POD basis for snapshots ysnap = G(usnap)
            % --------------------------------------------------------------------
            mode_d=0; %Snapshots for whole PDE
            timestart_podbasis = cputime;
            
            [ysnap,POD] = PODbasis(n,m,r,deltat,c,fxyt,au,usnap,M,K,Mt,y0,0,mode_d);
            time_podbasis = cputime - timestart_podbasis;
            
            % POD-DEIM basis for ysnap=ysnap^3
            %-------------------------------------------------------------------
            timestart_pod_deim_basis=cputime;
            mode_d=1; %Snapshots for non-linearity
 
            [~,U] = PODbasis(n,m, r_deim , deltat, 0 , 0 , 0 , 0 , M, 0 ,Mt, y0,ysnap, mode_d);
            
            [ P,p ] = DEIM( U, x1 , x2 , 2); 
            time_pod_deim_basis = cputime-timestart_pod_deim_basis;

            % Discretization of Omega by POD Galerkin method
            [Mpod,Kpod,y0pod,yqpod] = assem_POD(M,K,y0,yq,POD);

            % optimization by CG
            TOL2=2.83*1e-3;
            timestart_podopt = cputime;
            [yopt_pod_deim,uiopt_pod_deim,popt_pod_deim] = cg_pod_deim(n,r,r_deim, x1x2,tri,deltat,lambda,c,fxyt,ay,au,ua,ub,ustart,Mpod,Kpod,Mt,M,y0pod,yqpod,POD,U,P,linesearch_mode,cg_mode,TOL1,TOL2,s0,reset);
            time_podopt = cputime - timestart_podopt;
            time_pod_deim = cputime - timestart;
            disp(' ');
            disp(['CPU time - POD-DEIM optimization (incl. POD basis): ' num2str(time_pod_deim) 's']);
            disp(['   --> CPU time - POD basis: ' num2str(time_podbasis) 's']);
            disp(['   --> CPU time - POD-DEIM basis: ' num2str(time_pod_deim_basis) 's']);
            disp(['   --> CPU time - POD with DEIM Optimization: ' num2str(time_podopt) 's']);
            
 
            % sum(i=1:4)au_i(x)upod_i(t)
            uopt_pod_deim = zeros(m,n);
            for i = 1:n
                for j = 1:m
                    for k = 1:numcontrols
                        uopt_pod_deim(j,i) = uopt_pod_deim(j,i) + au(j,k)*uiopt_pod_deim(i,k);
                    end
                end
            end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the full equations for yopt, popt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(mode == 0 || mode == 1)
    
    [yopt_fem,popt_fem] = solve_full_eqn('FEM',n,m,uiopt_fem,y0,c,fxyt,au,M,K,deltat,yq,ay,mode,pod_mode);
        
end

if(mode == 0 || mode == 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     if pod_mode==0 || pod_mode==1

        [yfopt_pod,pf_pod] = solve_full_eqn('POD',n,m,uiopt_pod,y0,c,fxyt,au,M,K,deltat,yq,ay,mode,pod_mode);
        
     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POD-DEIM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     if pod_mode==0 || pod_mode==2

         [yfopt_pod_deim,pf_pod_deim] = solve_full_eqn('POD-DEIM',n,m,uiopt_pod_deim,y0,c,fxyt,au,M,K,deltat,yq,ay,mode,pod_mode);

     end
         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation Time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(' ');
    if(mode == 1 || mode == 0)
    disp(['CPU time - FEM optimization: ' num2str(time_fem) 's']);
    end
    if((mode == 2 || mode == 0) && (pod_mode==1 || pod_mode==0))
    disp(['CPU time - POD optimization (incl. POD basis): ' num2str(time_pod) 's']);
    end
    if((mode == 2 || mode == 0) && (pod_mode==2 || pod_mode==0))
    disp(['CPU time - POD-DEIM optimization (incl. POD basis): ' num2str(time_pod_deim) 's']);
    end
    disp(' ');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    minvalue = objvalue(uiopt,yopt,yq,ay,lambda,M,Mt);
    disp(' ');
    disp(['Objective function value: ' num2str(minvalue)]);

if(mode == 1 || mode == 0)
    minvalue_fem = objvalue(uiopt_fem,yopt_fem,yq,ay,lambda,M,Mt);
    disp(' ');
    disp(['Objective function value - FEM: ' num2str(minvalue_fem)]);
end
if((mode == 2 || mode == 0) && (pod_mode==1 || pod_mode==0))
    minvalue_pod = objvalue(uiopt_pod,yopt_pod,yq,ay,lambda,M,Mt);
    disp(' ');
    disp(['Objective function value - POD: ' num2str(minvalue_pod)]);
end
if((mode == 2 || mode == 0) && (pod_mode==2 || pod_mode==0))
    minvalue_pod_deim = objvalue(uiopt_pod_deim,yopt_pod_deim,yq,ay,lambda,M,Mt);
    disp(' ');
    disp(['Objective function value - POD-DEIM: ' num2str(minvalue_pod_deim)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(mode == 0 || mode == 1)
    
    disp_numerical_error('FEM',numcontrols,uiopt,uiopt_fem,yopt,yopt_fem,popt,popt_fem,M,deltat)  
        
end
% POD 
if((mode == 2 || mode == 0) && (pod_mode==1 || pod_mode==0))
    
    disp_numerical_error('POD',numcontrols,uiopt,uiopt_pod,yopt,yopt_pod,popt,popt_pod,M,deltat)  
    
end
% POD-DEIM
if((mode == 2 || mode == 0) && (pod_mode==2 || pod_mode==0))
        
  disp_numerical_error('POD-DEIM',numcontrols,uiopt,uiopt_pod_deim,yopt,yopt_pod_deim,popt,popt_pod_deim,M,deltat)  
      
end
% POD vs. FEM
if(mode == 0 && (pod_mode == 1 || pod_mode==0))
    
    disp_num_distance('POD','FEM',numcontrols, uiopt_fem,uiopt_pod,yopt_fem,yopt_pod,popt_fem,popt_pod,deltat,M)
    
end
%POD-DEIM vs. FEM
if(mode == 0 && (pod_mode == 2 || pod_mode==0))
    
    disp_num_distance('POD-DEIM','FEM',numcontrols, uiopt_fem,uiopt_pod_deim,yopt_fem,yopt_pod_deim,popt_fem,popt_pod_deim,deltat,M)
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical a-posteriori estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(mode == 1 && errormode == 1)
    if(c == 0)
        timestart = cputime;
        errana = error_analytical_linear(uiopt_fem, popt_fem,ua,ub,lambda,Mt,V1,V2,V3,V4,v);
        time_errana = cputime - timestart;            
        disp(' ');
        disp(['CPU time - analytical error estimation (linear): ' num2str(time_errana) 's']);
        disp(' ');
        disp(['Analytical error estimate of uopt_fem (linear) - FEM: ' num2str(errana)]);
    else
        timestart = cputime;
        [errana,mineig,cpu_hess] = error_analytical_nonlinear(uiopt_fem,yopt_fem,popt_fem,ua,ub,c,fxyt,y0,au,M,K,lambda,deltat,V1,V2,V3,V4,v,Mt,x1x2,tri);
        time_errana = cputime - timestart;            
        disp(' ');
        disp(['CPU time - complete analytical error estimation: ' num2str(time_errana) 's']);
        disp(['   --> CPU time - Establishing Reduced Hession: ' num2str(cpu_hess) 's']);
        disp(' ');
        disp(['Analytical error estimate of uopt_fem - FEM: ' num2str(errana)]);
        disp(' ');
        disp(['   --> smallest eigenvalue of reduced Hessian: ' num2str(mineig)]);
    end
end


if((mode == 2 || mode == 0) && errormode == 1)
    if(c == 0)
        timestart = cputime;
        errana = error_analytical_linear(uiopt_pod, popt_pod,ua,ub,lambda,Mt,V1,V2,V3,V4,v);
        time_errana = cputime - timestart;            
        disp(' ');
        disp(['CPU time - analytical error estimation (linear): ' num2str(time_errana) 's']);
        disp(' ');
        disp(['Analytical error estimate of uopt_pod (linear) - POD: ' num2str(errana)]);
    else
        timestart = cputime;
        [errana,mineig,cpu_hess] = error_analytical_nonlinear(uiopt_pod,yopt_pod,popt_pod,ua,ub,c,fxyt,y0,au,M,K,lambda,deltat,V1,V2,V3,V4,v,Mt,x1x2,tri);
        time_errana = cputime - timestart;            
        disp(' ');
        disp(['CPU time - complete analytical error estimation: ' num2str(time_errana) 's']);
        disp(['   --> CPU time - Establishing Reduced Hession: ' num2str(cpu_hess) 's']);
        disp(' ');
        disp(['Analytical error estimate of uopt_pod - POD: ' num2str(errana)]);
        disp(' ');
        disp(['   --> smallest eigenvalue of reduced Hessian: ' num2str(mineig)]);
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphical Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% %     
% %  Optimal solutions
% % ------------------------------------------------------------------------
% %     
%     % optimal control u1
%     figure;
%     subplot(2,2,1);
%     plot(t,uiopt(:,1),'k','Linewidth',2);
%     xlabel('time axis')
%     title('Optimal Control u_{opt}^1')
%     legend('u_1^{opt}',0);
%     
% 
%     % optimal control u2
%     subplot(2,2,2);
%     plot(t,uiopt(:,2),'k','Linewidth',2);
%     xlabel('time axis')
%     title('Optimal Control u_{opt}^2')
%     legend('u_2^{opt}',0);
%    
% 
%     % optimal control u3
%     subplot(2,2,3);
%     plot(t,uiopt(:,3),'k','Linewidth',2);
%     xlabel('time axis')
%     title('Optimal Control u_{opt}^3')
%     legend('u_3^{opt}',0);
%     
% 
%     % optimal control u4
%     subplot(2,2,4);
%     plot(t,uiopt(:,4),'k','Linewidth',2);
%     xlabel('time axis')
%     title('Optimal Control u_{opt}^4')
%     legend('u_4^{opt}',0);
%     
%     % time constant optimal state
%     figure;
%     subplot(1,2,1);
%     pdesurf(x1x2,tri,yopt(:,1));
%     xlabel('x1 axis');
%     ylabel('x2 axis');
%     title('Optimal State y_{opt} (time constant)');
%     
% 
%     % optimal adjoint state p
%     subplot(1,2,2);
%     pdesurf(x1x2,tri,popt(:,floor(n/2)));
%     xlabel('x1 axis');
%     ylabel('x2 axis');
%     title('Optimal Adjoint State p_{opt} at time T/2');
 % Comparisons optimal vs. FEM
% ------------------------------------------------------------------------
%   
if ( mode == 1 || mode == 0)
    % optimal control uiopt vs. uiopt_fem
    figure;
    for i = 1:numcontrols
        subplot(2,2,i);
        plot(t,uiopt(:,i),'k',t,uiopt_fem(:,i),'r','Linewidth',2);
        xlabel('time axis')
        title(['Optimal Control u_{opt}^' num2str(i) ' vs. u_{fem}^' num2str(i)])
        legend(['u_{opt}^' num2str(i)],['u_{fem}^' num2str(i)],0);
         
    end
    
    % optimality systems uifemi <---> pfem
  
    

     optimality_test('FEM',m,n,t,V1,V2,V3,V4,popt_fem,uiopt_fem)
    
end


% Comparisons optimal vs. POD
% ------------------------------------------------------------------------
%   
if ( (mode == 2 || mode == 0) && (pod_mode==1 || pod_mode==0))
    % optimal control uiopt vs. uiopt_pod
    figure;
    for i = 1:numcontrols
        subplot(2,2,i);
        plot(t,uiopt(:,i),'k',t,uiopt_pod(:,i),'c','Linewidth',2);
        xlabel('time axis')
        title(['Optimal Control u_{opt}^' num2str(i) ' vs. u_{pod}^' num2str(i)])
        legend(['u_{opt}^' num2str(i)],['u_{pod}^' num2str(i)],0);
        
    end
    
    % with u compute by POD
    % optimality systems uipodi <---> ppod

    optimality_test('POD',m,n,t,V1,V2,V3,V4,popt_pod,uiopt_pod)
   
    
    % with u compute by POD and y,p by the full equation
    % optimality systems uipodi <---> ppod
    
    optimality_test('PODfull',m,n,t,V1,V2,V3,V4,pf_pod,uiopt_pod)
   
end

% Comparisons optimal vs. POD-DEIM
% ------------------------------------------------------------------------
%   
if ( (mode == 2 || mode == 0) && (pod_mode==2 || pod_mode==0))
    % optimal control uiopt vs. uiopt_pod_deim
    figure;
    for i = 1:numcontrols
        subplot(2,2,i);
        plot(t,uiopt(:,i),'k',t,uiopt_pod_deim(:,i),'c','Linewidth',2);
        xlabel('time axis')
        title(['Optimal Control u_{opt}^' num2str(i) ' vs. u_{pod-deim}^' num2str(i)])
        legend(['u_{opt}^' num2str(i)],['u_{pod-deim}^' num2str(i)],0);
        
    end
    
    % optimality systems uipodi <---> ppod
    
    optimality_test('POD-DEIM', m,n,t,V1,V2,V3,V4,popt_pod_deim,uiopt_pod_deim)

    
    % with u compute by POD and y,p by the full equation
    % optimality systems uipodi <---> ppod
    
    optimality_test('POD-DEIMfull',m,n,t,V1,V2,V3,V4,pf_pod_deim,uiopt_pod_deim)
    
    
end

% Comparisons FEM vs. POD
% ------------------------------------------------------------------------
%   
% if mode == 0
%     % uiopt_fem vs. uiopt_pod
%     for i = 1:numcontrols
%         figure;
%         plot(t,uiopt_fem(:,i),'r',t,uiopt_pod(:,i),'c','Linewidth',2);
%         xlabel('time axis')
%         title(['u_{fem}^' num2str(i) ' vs. u_{pod}^' num2str(i)])
%         legend(['u_{fem}^' num2str(i)],['u_{pod}^' num2str(i)],0);
%         grid on; 
%     end
% end

% weight functions au_i
figure;
pdesurf(x1x2,tri,au(:,1)+au(:,2)+au(:,3)+au(:,4));
title('Weight functions au_i(x)');
xlabel('x1 axis');
ylabel('x2 axis');
axis([0 pi 0 pi 0 10]);
grid on;


disp(' ');
disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
disp(' ');
disp(' ');
     
