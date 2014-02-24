function example2D_ssn(mmin,n,r,r_deim,mode,errormode,pod_mode)

% =========================================================================
% Author: Sophia Kohle, Technische Universität Berlin
% =========================================================================
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Task: Computation of an optimal control for the problem
%
% min 1/2 ||y-yq||^2 + lambda/2 sum(i=1:4) ||u_i||^2 + (ay,y(T))
% s.t.
%       y_t(x,t) - /\y(x,t) + c y^3(x,t) + d(x,t) = sum(i=1:4)au_i(x) u_i(t) 
%       y(x,0) = y_0(x) 
%       dy/dnu(x,t) = 0
%
%       ua <= u_i(t) <= ub     , i = 1:4
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% 
% Used algorithms:  Semi smooth Newton method, primal-dual active-set method, finite
%                   element discretization in space, semi-implicit Euler
%                   method in time
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Input parameters:
%              
% mmin       % minimal number of finite elements in space discretization
% n          % number of time instances in [0,T]
% r          % rank of POD basis
% mode       % 0 = Computation by FE method and by POD method
%            % 1 = Computation by FE method
%            % 2 = Computation by POD method
% errormode  % 0 = no computation of analytical a-posteriori-estimate
%            % 1 = computation of analytical a-posteriori-estimate
%                
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
if(mode == 0 || mode == 1)
    
    % optimization by SQP
        disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
        disp('    Optimization by using FEM');
        disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
        timestart = cputime;
    [yopt_fem,uiopt_fem,popt_fem] = ssn_fem(n,m,x1x2,tri,deltat,lambda,c,fxyt,ay,au,ua,ub,ustart,ystart,pstart,M,K,y0,yq,V1,V2,V3,V4,v,Mt);
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

    % Discretization of (0,L) by POD Galerkin method
    [Mpod,Kpod,y0pod,yqpod] = assem_POD(M,K,y0,yq,POD);
    
    % optimization by SQP
        
        timestart_podopt = cputime;
        [yopt_pod,uiopt_pod,popt_pod] = ssn_pod(n,r,x1x2,tri,deltat,lambda,c,POD'*M*fxyt,POD'*M*ay,POD'*M*au,ua,ub,ustart,POD'*M*ystart,POD'*M*pstart,Mpod,Kpod,Mt,y0pod,yqpod,POD);
        time_podopt = cputime - timestart_podopt;
        time_pod = cputime - timestart;
        disp(' ');
        disp(['CPU time - POD optimization (incl. POD basis): ' num2str(time_pod) 's']);
        disp(['   --> CPU time - POD basis: ' num2str(time_podbasis) 's']);
        disp(['   --> CPU time - POD Optimization: ' num2str(time_podopt) 's']);

        % sum(i=1:4)au_i(x)upod_i(t)
        uopt_pod = zeros(m,n);
        for i = 1:n
            for j = 1:m
                for k = 1:numcontrols
                    uopt_pod(j,i) = uopt_pod(j,i) + au(j,k)*uiopt_pod(i,k);
                end
            end
        end
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
           
            [ P,p ] = DEIM( U, x1 , x2 , 2); %***
            time_pod_deim_basis = cputime-timestart_pod_deim_basis;

            % Discretization of Omega by POD Galerkin method
            [Mpod,Kpod,y0pod,yqpod] = assem_POD(M,K,y0,yq,POD);

            % optimization by steepest descent
            TOL2 = 5*1e-5;
            timestart_podopt = cputime;
            [yopt_pod_deim,uiopt_pod_deim,popt_pod_deim] = ssn_pod_deim(n,r,x1x2,tri,deltat,lambda,c,POD'*M*fxyt,POD'*M*ay,POD'*M*au,ua,ub,ustart,POD'*M*ystart,POD'*M*pstart,Mpod,Kpod,Mt,y0pod,yqpod,POD,U,P,M);
            time_podopt = cputime - timestart_podopt;
            time_pod_deim = cputime - timestart;
            disp(' ');
            disp(['CPU time - POD-DEIM optimization (incl. POD basis): ' num2str(time_pod_deim) 's']);
            disp(['   --> CPU time - POD basis: ' num2str(time_podbasis) 's']);
            disp(['   --> CPU time - POD-DEIM basis: ' num2str(time_pod_deim_basis) 's']);
            disp(['   --> CPU time - POD with DEIM Optimization: ' num2str(time_podopt) 's']);


            % sum(i=1:4)au_i(x)upoddeim_i(t)
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphical Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
     

save example2_results.mat;
