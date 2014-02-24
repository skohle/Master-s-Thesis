function disp_num_distance(string1,string2,numcontrols,uiopt1,uiopt2,yopt1,yopt2,popt1,popt2,deltat,M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displays the numerical error of u, y and p of string1 and string2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:numcontrols
        % ui_opt vs. ui_pod
        errnum(i) = error_numerical_0T(uiopt1(:,i), uiopt2(:,i), deltat);
    end
    % ||u||_L²(0,T;R4),   u = (u1,...,u4)'
    errnum = sqrt(sum(errnum.*errnum));
    disp(' ');
    disp(['Numerical distance between ' string1 ' and ' string2 ' of u: ' num2str(errnum)]);
    
    % ||y||_L²(Q)
    errnum = error_numerical_Q(yopt1, yopt2, M, deltat);
    disp(['                           ' string1 ' and ' string2 ' of y: ' num2str(errnum)]);
    % ||p||_L²(Q)
    errnum = error_numerical_Q(popt1, popt2, M, deltat);
    disp(['                           ' string1 ' and ' string2 ' of p: ' num2str(errnum)]);
    disp('   ');