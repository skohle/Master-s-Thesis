function example2D_ml_BFGS(r,r_deim,mode,errormode,cg_mode,linesearch_mode,pod_mode,mminvec,nvec,bfgs_mode)

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
% 
%
% =========================================================================
% 
%
%
% Parameters
%----------------------------------------------------------------------

T = 1;                              % Final time point
numcontrols = 4;             % number of control functions

 
lambda = 0.01;               % Tikhonov parameter (small)
c = 1;                       % parameter c in PDE

TOL1 = 5*1e-5;
TOL2=3.6*1e-3;
TOL=TOL2;
s0=1;                       % Start stepsize
reset=5;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Computation by FE method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(mode == 0 || mode == 1)

    time_fem_start=cputime;
    TOL1=5*1e-3;

    % First solving                        
    ustart = sparse(nvec(1),numcontrols);
    mmin=mminvec(1);
    n=nvec(1);
    deltat = T/(n-1);
    t = (0:deltat:T)';
    disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
    disp(['    Multilevel-Optimization 1 ']);
    disp(['    Spatial discretisation mmin=' num2str(mmin)]);
    disp(['    Time discretisation n=' num2str(n)]);
    disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
    
    [yopt_fem,uiopt_fem,popt_fem,m,x1,x2,au,M,V1,V2,V3,V4,v,yq,Mt,ay] = multilevel_bfgs_fem(mmin,n,  mode, linesearch_mode, ustart , t, lambda, c, deltat,T,TOL1,TOL2,s0,bfgs_mode); 
    
    % plot the projection
    % --------------------------------------------------------------
    % optimality systems uifemi <---> pfem
    optimality_test('FEM',m,n,t,V1,V2,V3,V4,popt_fem,uiopt_fem)
    
    
    for i=2:length(mminvec);
 
    % Interpolation on the new time-grid
           t_old = t;
           mmin=mminvec(i);
           n=nvec(i);
           deltat = T/(n-1);
           t = (0:deltat:T)';
           TOL1 = 5*1e-5;
           disp('                                                     ');
           disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
           disp(['    Multilevel-Optimization ' num2str(i) ]);
           disp(['    Spatial discretisation mmin=' num2str(mmin)]);
           disp(['    Time discretisation n=' num2str(n)]);
           disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
           disp('                                                     ')
           ustart = sparse(n,numcontrols);
           ustart(:,1)=interp1(t_old,uiopt_fem(:,1),t);
           ustart(:,2)=interp1(t_old,uiopt_fem(:,2),t);
           ustart(:,3)=interp1(t_old,uiopt_fem(:,3),t);
           ustart(:,4)=interp1(t_old,uiopt_fem(:,4),t);
           [yopt_fem,uiopt_fem,popt_fem,m,x1,x2,au,M,V1,V2,V3,V4,v,yq,Mt,ay,y0,fxyt,K,tri] =multilevel_bfgs_fem(mmin, n,mode, linesearch_mode, ustart, t, lambda, c, deltat,T,TOL1,TOL2,s0,bfgs_mode); 
    %    end
    

    % plot the projection
    % --------------------------------------------------------------
    % optimality systems uifemi <---> pfem
    optimality_test('FEM',m,n,t,V1,V2,V3,V4,popt_fem,uiopt_fem)
    
    end
time_fem_end=cputime-time_fem_start;

disp(['CPU time - Multilevel FEM optimization: ' num2str(time_fem_end) 's']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2a) Computation by POD method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(mode == 0 || mode == 2)
    
    if pod_mode==0 || pod_mode==1
        

        time_pod_start=cputime;
        TOL1=5*1e-5;
        TOL2=3.23*1e-3;
        % First solving                        
        ustart = sparse(nvec(1),numcontrols);
        mmin=mminvec(1);
        n=nvec(1);
        deltat = T/(n-1);
        t = (0:deltat:T)';
        usnap(:,1) = 0.5*ones(n,1);  % controls for generation of snapshots
        usnap(:,2) = -0.5*ones(n,1);
        usnap(:,3) = 0.5*ones(n,1);
        usnap(:,4) = -0.5*ones(n,1);
        
        disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
        disp(['    Multilevel-Optimization 1 by using POD with basis rank r = ' num2str(r)]);
        disp(['    Spatial discretisation mmin=' num2str(mmin)]);
        disp(['    Time discretisation n=' num2str(n)]);
        disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
        

        [yopt_pod,uiopt_pod,popt_pod,m,x1,x2,au,M,V1,V2,V3,V4,v,yq,Mt,ay,y0,fxyt,K,tri] = multilevel_bfgs_pod(mmin, n, r,mode,linesearch_mode,pod_mode, ustart , t, lambda, c, deltat,T,TOL1,TOL2,s0,usnap,1,bfgs_mode); 
        
        % plot the projection
            % --------------------------------------------------------------
            % optimality systems uifemi <---> pfem
            optimality_test('POD',m,n,t,V1,V2,V3,V4,popt_pod,uiopt_pod)

        for i=2:length(mminvec)
   
            % Interpolation on the new time-grid
                   t_old = t;
                   mmin=mminvec(i);
                   n=nvec(i);
                   deltat = T/(n-1);
                   t = (0:deltat:T)';
                   disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
                   disp(['    Multilevel-Optimization '  num2str(i) ' by using POD with basis rank r = ' num2str(r)]);
                   disp(['    Spatial discretisation mmin=' num2str(mmin)]);
                   disp(['    Time discretisation n=' num2str(n)]);
                   disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
                    s0=4;
                   ustart = sparse(n,numcontrols);
                   ustart(:,1)=interp1(t_old,uiopt_pod(:,1),t);
                   ustart(:,2)=interp1(t_old,uiopt_pod(:,2),t);
                   ustart(:,3)=interp1(t_old,uiopt_pod(:,3),t);
                   ustart(:,4)=interp1(t_old,uiopt_pod(:,4),t);
                   usnap=ustart;
                   [yopt_pod,uiopt_pod,popt_pod,m,x1,x2,au,M,V1,V2,V3,V4,v,yq,Mt,ay,y0,fxyt,K,tri] =multilevel_bfgs_pod(mmin, n, r, mode, linesearch_mode,pod_mode, ustart, t, lambda, c, deltat,T,TOL1,TOL2,s0,usnap,i,bfgs_mode); 
            %    end


            % plot the projection
            % --------------------------------------------------------------
            % optimality systems uifemi <---> pfem
            optimality_test('POD',m,n,t,V1,V2,V3,V4,popt_pod,uiopt_pod)

        end
        time_pod_end=cputime-time_pod_start;

        disp(['CPU time - multilevel - POD optimization: ' num2str(time_pod_end) 's']);

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2b) Computation by POD-DEIM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if pod_mode==0 || pod_mode==2
            time_pod_deim_start=cputime;
            clear usnap
            % First solving                        
            ustart = sparse(nvec(1),numcontrols);
   
            mmin=mminvec(1);
            n=nvec(1);
            deltat = T/(n-1);
            t = (0:deltat:T)';
            usnap(:,1) = 0.5*ones(n,1);  % controls for generation of snapshots
            usnap(:,2) = -0.5*ones(n,1);
            usnap(:,3) = 0.5*ones(n,1);
            usnap(:,4) = -0.5*ones(n,1);

            disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
            disp(['    Multilevel-Optimization 1 by using POD-DEIM with basis rank r = ' num2str(r)]);
            disp(['    Spatial discretisation mmin=' num2str(mmin)]);
            disp(['    Time discretisation n=' num2str(n)]);
            disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
            TOL1=1*1e-5;
            TOL2=2*1e-3;

            [yopt_pod_deim,uiopt_pod_deim,popt_pod_deim,m,x1,x2,au,M,V1,V2,V3,V4,v,yq,Mt,ay,y0,fxyt,K,tri] = multilevel_bfgs_pod_deim(mmin, n, r,r_deim, mode, linesearch_mode,pod_mode, ustart ,  lambda, c, deltat,TOL1,TOL2,s0,usnap,bfgs_mode); 
            
            
            % plot the projection
                % --------------------------------------------------------------
                % optimality systems uifemi <---> pfem
                optimality_test('POD-DEIM', m,n,t,V1,V2,V3,V4,popt_pod_deim,uiopt_pod_deim)
                
                
            for i=2:length(mminvec)
         
                % Interpolation on the new time-grid
                if i==2
                    TOL2=5*1e-3;
                end
                if i==3
                    TOL2=3.5*1e-3;
                    %s0=2;
                end
                       t_old = t;
                       mmin=mminvec(i);
                       n=nvec(i);
                       deltat = T/(n-1);
                       t = (0:deltat:T)';
                       disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');
                       disp(['    Multilevel-Optimization '  num2str(i) ' by using POD-DEIM with basis rank r = ' num2str(r)]);
                       disp(['    Spatial discretisation mmin=' num2str(mmin)]);
                       disp(['    Time discretisation n=' num2str(n)]);
                       disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');

                       ustart = sparse(n,numcontrols);
                       ustart(:,1)=interp1(t_old,uiopt_pod_deim(:,1),t);
                       ustart(:,2)=interp1(t_old,uiopt_pod_deim(:,2),t);
                       ustart(:,3)=interp1(t_old,uiopt_pod_deim(:,3),t);
                       ustart(:,4)=interp1(t_old,uiopt_pod_deim(:,4),t);
                       
                       usnap=ustart;
                       
                       [yopt_pod_deim,uiopt_pod_deim,popt_pod_deim,m,x1,x2,au,M,V1,V2,V3,V4,v,yq,Mt,ay,y0,fxyt,K,tri] =multilevel_bfgs_pod_deim(mmin, n, r,r_deim, mode, linesearch_mode,pod_mode, ustart, lambda, c, deltat,TOL1,TOL2,s0,usnap,bfgs_mode); 
                %    end


                % plot the projection
                % --------------------------------------------------------------
                % optimality systems uifemi <---> pfem
                optimality_test('POD-DEIM', m,n,t,V1,V2,V3,V4,popt_pod_deim,uiopt_pod_deim)
            end
            time_pod_deim_end=cputime-time_pod_deim_start;

            disp(['CPU time - multilevel - POD-DEIM optimization: ' num2str(time_pod_deim_end) 's']);


        end
end



ua = -ones(n,numcontrols);   % lower control bound
ub = ones(n,numcontrols);    % upper control bound 


% constructed optimal solutions
% ------------------------------------------------------------------------
 
% optimal state y
yopt = zeros(m,n);
for i = 1:n
    yopt(:,i) = cos(x1).*cos(x2);
end
% optimal adjoint state
popt = zeros(m,n);
for i = 1:n
    for k = 1:m
        popt(k,i) = lambda*t(i)^2*cos(x1(k))*cos(x2(k));
    end
end
% optimal control
    pxy = cos(x1).*cos(x2);
    Aaup = au'*M*pxy;
uiopt = zeros(n,numcontrols);
for i = 1:numcontrols
    uiopt(:,i) = min(ub(:,i),max(ua(:,i),-t.^2*Aaup(i)));
end

    uopt = zeros(m,n);
    for i = 1:n
        for j = 1:m
            for k = 1:numcontrols
                uopt(j,i) = uopt(j,i) + au(j,k)*uiopt(i,k);
            end
        end
    end
x1x2(1,:) = x1';
x1x2(2,:) = x2';
   
   
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
    if ( mode == 1 || mode == 0)
    disp(['CPU time - multilevel - FEM optimization: ' num2str(time_fem_end) 's']);
    end
    if ( (mode == 2 || mode == 0) && (pod_mode==1 || pod_mode==0))
    disp(['CPU time - multilevel -  POD optimization (incl. POD basis): ' num2str(time_pod_end) 's']);
    end
    if ( (mode == 2 || mode == 0) && (pod_mode==2 || pod_mode==0))
    disp(['CPU time - multilevel - POD-DEIM optimization (incl. POD basis): ' num2str(time_pod_deim_end) 's']);
    end
    disp(' ');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minvalue = objvalue(uiopt,yopt,yq,ay,lambda,M,Mt);
% disp(' ');
% disp(['Minimal objective function value: ' num2str(minvalue)]);

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
%  Optimal solutions
% ------------------------------------------------------------------------
%     
    % optimal control u1
    figure;
    subplot(2,2,1);
    plot(t,uiopt(:,1),'k','Linewidth',2);
    xlabel('time axis')
    title('Optimal Control u_{opt}^1')
    legend('u_1^{opt}',0);
    

    % optimal control u2
    subplot(2,2,2);
    plot(t,uiopt(:,2),'k','Linewidth',2);
    xlabel('time axis')
    title('Optimal Control u_{opt}^2')
    legend('u_2^{opt}',0);
   

    % optimal control u3
    subplot(2,2,3);
    plot(t,uiopt(:,3),'k','Linewidth',2);
    xlabel('time axis')
    title('Optimal Control u_{opt}^3')
    legend('u_3^{opt}',0);
    

    % optimal control u4
    subplot(2,2,4);
    plot(t,uiopt(:,4),'k','Linewidth',2);
    xlabel('time axis')
    title('Optimal Control u_{opt}^4')
    legend('u_4^{opt}',0);
    
    % time constant optimal state
    figure;
    subplot(1,2,1);
    pdesurf(x1x2,tri,yopt(:,1));
    xlabel('x1 axis');
    ylabel('x2 axis');
    title('Optimal State y_{opt} (time constant)');
    

    % optimal adjoint state p
    subplot(1,2,2);
    pdesurf(x1x2,tri,popt(:,floor(n/2)));
    xlabel('x1 axis');
    ylabel('x2 axis');
    title('Optimal Adjoint State p_{opt} at time T/2');
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
x1x2(1,:) = x1';
x1x2(2,:) = x2';
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