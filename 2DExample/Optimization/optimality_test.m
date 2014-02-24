function optimality_test(string,m,n,t,V1,V2,V3,V4,p,u)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  plot of the optimality test for given u and p
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    popt_vec = reshape(p,m*n,1);
    figure;
    subplot(2,2,1);
    plot(t,-V1*popt_vec,'--',t,u(:,1),'k','Linewidth',2);
    title(strcat(string , ' Optimality test u^1'));
    xlabel('time axis')
    legend('-1/lambda V_1^*p','u^1',0)
    subplot(2,2,2);
    plot(t,-V2*popt_vec,'--',t,u(:,2),'k','Linewidth',2);
    title(strcat(string , ' Optimality test u^2'));
    xlabel('time axis')
    legend('-1/lambda V_2^*p','u^2',0)
    subplot(2,2,3);
    plot(t,-V3*popt_vec,'--',t,u(:,3),'k','Linewidth',2);
    title(strcat(string ,  ' Optimality test u^3'));
    xlabel('time axis')
    legend('-1/lambda V_3^*p','u^3',0)
    subplot(2,2,4);
    plot(t,-V4*popt_vec,'--',t,u(:,4),'k','Linewidth',2);
    title(strcat(string ,  ' Optimality test u^4'));
    xlabel('time axis')
    legend('-1/lambda V_4^*p','u^4',0)

