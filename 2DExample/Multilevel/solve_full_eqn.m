function [yopt,popt]=solve_full_eqn(string,n,m,ui,y0,c,fxyt,au,M,K,deltat,yq,ay,mode,pod_mode)


timestart = cputime;
    yopt = state_equation(n,m,ui,y0,c,fxyt,au,M,K,deltat);
        time_y = cputime - timestart;
        disp(' ');
        disp(['CPU time - Computation of full ' string '-y: ' num2str(time_y) 's']);
    
        timestart = cputime;
    popt = adjoint_equation(n,m,yopt,deltat,yq,c,ay,M,K);
        time_p = cputime - timestart;
        disp(['CPU time - Computation of full ' string '-p: ' num2str(time_p) 's']);
    