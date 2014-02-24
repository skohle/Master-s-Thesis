function disp_numerical_error(string,numcontrols,uiopt,ui,yopt,y,popt,p,M,deltat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Displays the numerical error of u, p, y to the exact solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:numcontrols
        % ui_opt vs. ui
        errnum(i) = error_numerical_0T(uiopt(:,i), ui(:,i), deltat);
end
    % ||u||_L²(0,T;R4),   u = (u1,...,u4)'
    errnum = sqrt(sum(errnum.*errnum));
    disp(' ');
    disp(['Numerical error ' string '(ui to uiopt): ' num2str(errnum)]);
    
    % ||y||_L²(Q)
    errnum = error_numerical_Q(yopt, y, M, deltat);
    disp(['                ' string '(y  to yopt): ' num2str(errnum)]);
    % ||p||_L²(Q)
    errnum = error_numerical_Q(popt, p, M, deltat);
    disp(['                ' string ' (p  to popt): ' num2str(errnum)]);
    disp('   ');
    