function [yoptvec, uopt, poptvec] = active_set(ua,ub,G1,G21,G22,G23,G24,g,ADG1,ADG2,adg,V1,V2,V3,V4,v,m,n,lambda,ystart,pstart)

% =========================================================================
% Author: Sophia Kohle, Technische Universität Berlin
% =========================================================================
%
% Task: Solution of a linear-quadratic optimal control problem with control
%       constraints
%
%       ua <= u(t) <= ub
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Numerical methods:    One iteration of primal-dual Active Set Method
%              
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Input parameters:
%
% m                         % number of finite elements in [0,L]
% n                         % number of time instances in [0,T]
% ua                        % lower control bound
% ub                        % upper control bound
% G1,G21,G22,G23,G24,g      % matrices for state equation
% ADG1, ADG2, adg           % matrices for adjoint equation
% V1,V2,V3,V4,v             % matrices for projection formula
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Output:   optimal solution set (uopt, yopt,popt)
% =========================================================================

% Parameters

numcontrols = 4;
cmue = lambda;

% Step (0)
% ------------------------------------------------------------------------

p = pstart;
y = ystart;
u = zeros(n,numcontrols);
u(:,1) = -V1*p-v;
u(:,2) = -V2*p-v;
u(:,3) = -V3*p-v;
u(:,4) = -V4*p-v;

mue = cmue*sparse(n,numcontrols);
    
    % Step (1) - Computing of active/inactive sets
    % ------------------------------------------------------------------------
    Aaset = [];
    Abset = [];
    Iset = [];

    for k = 1:numcontrols
        Aasetk = [];
        Absetk = [];
        Isetk = [];
        for i = 1:n
            % Active Set Ak_a
            if u(i,k) - ua(i,k) + mue(i,k) <= 0
                Aasetk = [Aasetk,i];
            % Active Set Ak_b
            elseif u(i,k) - ub(i,k) + mue(i,k) >= 0
                Absetk = [Absetk,i];
            else
            % Inactive Set Ik
                Isetk = [Isetk,i];
            end
        end

        if k == 1
            Aaset(:,:,k) = Aasetk;
            Abset(:,:,k) = Absetk;
            Iset(:,:,k) = Isetk;
        else
            % Aa
            diffsizeAa = size(Aasetk,2) - size(Aaset(:,:,k-1),2);
            if diffsizeAa == 0
                Aaset(:,:,k) = Aasetk;
            elseif diffsizeAa < 0
                Aaset(:,:,k) = [Aasetk,zeros(1,-diffsizeAa)];
            elseif diffsizeAa > 0
                Aasethelp = Aaset;
                Aaset = [];
                for l = 1:k-1
                    Aaset(:,:,l) = [Aasethelp(:,:,l),zeros(1,diffsizeAa)];
                end
                Aaset(:,:,k) = Aasetk;
            end
            % Ab
            diffsizeAb = size(Absetk,2) - size(Abset(:,:,k-1),2);
            if diffsizeAb == 0
                Abset(:,:,k) = Absetk;
            elseif diffsizeAb < 0
                Abset(:,:,k) = [Absetk,zeros(1,-diffsizeAb)];
            elseif diffsizeAb > 0
                Absethelp = Abset;
                Abset = [];
                for l = 1:k-1
                    Abset(:,:,l) = [Absethelp(:,:,l),zeros(1,diffsizeAb)];
                end
                Abset(:,:,k) = Absetk;
            end
            % I
            diffsizeI = size(Isetk,2) - size(Iset(:,:,k-1),2);
            if diffsizeI == 0
                Iset(:,:,k) = Isetk;
            elseif diffsizeI < 0
                Iset(:,:,k) = [Isetk,zeros(1,-diffsizeI)];
            elseif diffsizeI > 0
                Isethelp = Iset;
                Iset = [];
                for l = 1:k-1
                    Iset(:,:,l) = [Isethelp(:,:,l),zeros(1,diffsizeI)];
                end
                Iset(:,:,k) = Isetk;
            end

        end
    end



    % matrices for active sets
    Eaktivminus1 = sparse(n,n);
    Eaktivminus2 = sparse(n,n);
    Eaktivminus3 = sparse(n,n);
    Eaktivminus4 = sparse(n,n);
    Eaktivplus1 = sparse(n,n);
    Eaktivplus2 = sparse(n,n);
    Eaktivplus3 = sparse(n,n);
    Eaktivplus4 = sparse(n,n);
    Einaktiv1 = sparse(n,n);
    Einaktiv2 = sparse(n,n);
    Einaktiv3 = sparse(n,n);
    Einaktiv4 = sparse(n,n);

    
    % Eaktivminus1
    for i = 1:size(Aaset(:,:,1),2)
        if Aaset(1,i,1) == 0
            break;
        else
            Eaktivminus1(Aaset(1,i,1),Aaset(1,i,1))=1;
        end
    end
    % Eaktivplus1
    for i = 1:size(Abset(:,:,1),2)
        if Abset(1,i,1) == 0
            break;
        else
            Eaktivplus1(Abset(1,i,1),Abset(1,i,1))=1;
        end
    end    
    % Einaktiv1
    for i = 1:size(Iset(:,:,1),2)
        if Iset(1,i,1) == 0
            break;
        else
            Einaktiv1(Iset(1,i,1),Iset(1,i,1))=1;
        end
    end
    % Eaktiv1
    Eaktiv1 = Eaktivminus1+Eaktivplus1;
    % Eaktivminus2
    for i = 1:size(Aaset(:,:,2),2)
        if Aaset(1,i,2) == 0
            break;
        else
            Eaktivminus2(Aaset(1,i,2),Aaset(1,i,2))=1;
        end
    end
    % Eaktivplus2
    for i = 1:size(Abset(:,:,2),2)
        if Abset(1,i,2) == 0
            break;
        else
            Eaktivplus2(Abset(1,i,2),Abset(1,i,2))=1;
        end
    end    
    % Einaktiv2
    for i = 1:size(Iset(:,:,2),2)
        if Iset(1,i,2) == 0
            break;
        else
            Einaktiv2(Iset(1,i,2),Iset(1,i,2))=1;
        end
    end
    % Eaktiv2
    Eaktiv2 = Eaktivminus2+Eaktivplus2;
    % Eaktivminus3
    for i = 1:size(Aaset(:,:,3),2)
        if Aaset(1,i,3) == 0
            break;
        else
            Eaktivminus3(Aaset(1,i,3),Aaset(1,i,3))=1;
        end
    end
    % Eaktivplus3
    for i = 1:size(Abset(:,:,3),2)
        if Abset(1,i,3) == 0
            break;
        else
            Eaktivplus3(Abset(1,i,3),Abset(1,i,3))=1;
        end
    end    
    % Einaktiv3
    for i = 1:size(Iset(:,:,3),2)
        if Iset(1,i,3) == 0
            break;
        else
            Einaktiv3(Iset(1,i,3),Iset(1,i,3))=1;
        end
    end
    % Eaktiv3
    Eaktiv3 = Eaktivminus3+Eaktivplus3;
    % Eaktivminus4
    for i = 1:size(Aaset(:,:,4),2)
        if Aaset(1,i,4) == 0
            break;
        else
            Eaktivminus4(Aaset(1,i,4),Aaset(1,i,4))=1;
        end
    end
    % Eaktivplus4
    for i = 1:size(Abset(:,:,4),2)
        if Abset(1,i,4) == 0
            break;
        else
            Eaktivplus4(Abset(1,i,4),Abset(1,i,4))=1;
        end
    end    
    % Einaktiv4
    for i = 1:size(Iset(:,:,4),2)
        if Iset(1,i,4) == 0
            break;
        else
            Einaktiv4(Iset(1,i,4),Iset(1,i,4))=1;
        end
    end
    % Eaktiv4
    Eaktiv4 = Eaktivminus4+Eaktivplus4;
    
    
    
    % setting fixed variables of u
    uaktiv(:,1) = Eaktivminus1*ua(:,1) + Eaktivplus1*ub(:,1);
    uaktiv(:,2) = Eaktivminus2*ua(:,2) + Eaktivplus2*ub(:,2);
    uaktiv(:,3) = Eaktivminus3*ua(:,3) + Eaktivplus3*ub(:,3);
    uaktiv(:,4) = Eaktivminus4*ua(:,4) + Eaktivplus4*ub(:,4);
    
    % test
    if(~isequal(Eaktiv1+Einaktiv1,speye(n)))
        error('error initializing active sets');
    end
    if(~isequal(Eaktiv2+Einaktiv2,speye(n)))
        error('error initializing active sets');
    end
    if(~isequal(Eaktiv3+Einaktiv3,speye(n)))
        error('error initializing active sets');
    end
    if(~isequal(Eaktiv4+Einaktiv4,speye(n)))
        error('error initializing active sets');
    end
    
    
    
    % Step (2) - Solving LES for y, u, p
    % ------------------------------------------------------------------------
   
    % 1. rows for PDE constraints
    % G1 y = G2 u + g
    % G1 y = G2 (-V p - v) + g
    %                                               H(1:n*m,1:n*m) = G1;
    %                                               H(1:n*m,n*m+1:2*n*m) = G2*Einaktiv*V;
    %                                               h(1:n*m) = g + G2*uaktiv - G2*Einaktiv*v;
    G2EinaktivV = G21*Einaktiv1*V1+ G22*Einaktiv2*V2 + G23*Einaktiv3*V3 + G24*Einaktiv4*V4;
    G2uaktiv = G21*uaktiv(:,1) + G22*uaktiv(:,2) + G23*uaktiv(:,3) + G24*uaktiv(:,4);
    G2Einaktivv = (G21*Einaktiv1 + G22*Einaktiv2 + G23*Einaktiv3 + G24*Einaktiv4)*v;

    
    [iG1,jG1,sG1] = find(G1);
    [iG2,jG2,sG2] = find(G2EinaktivV);
    if ~isequal(jG2,[])
        jG2 = jG2 + n*m*ones(length(iG2),1);
    end
    % 2. rows for adjoint equation
    % ADG1 p = ADG2 y + adg
    %                                               H(n*m+1:2*n*m,1:n*m) = -ADG2;
    %                                               H(n*m+1:2*n*m,n*m+1:2*n*m) = ADG1;
    %                                               h(n*m+1:2*n*m) = adg;
    [iADG2,jADG2,sADG2] = find(-ADG2);
    if ~isequal(iADG2,[])
        iADG2 = iADG2 + n*m*ones(length(iADG2),1);
    end
    [iADG1,jADG1,sADG1] = find(ADG1);
    if ~isequal(iADG1,[])
        iADG1 = iADG1 + n*m*ones(length(iADG1),1);
        jADG1 = jADG1 + n*m*ones(length(jADG1),1);
    end

    % H, h
    i = [iG1;iG2;iADG2;iADG1];
    j = [jG1;jG2;jADG2;jADG1];
    s = [sG1;sG2;sADG2;sADG1];
    H = sparse(i,j,s, 2*n*m, 2*n*m);

    [ih1,jh1,sh1] = find(g + G2uaktiv - G2Einaktivv);
    [iadg,jadg,sadg] = find(adg);
    ih = [ih1;n*m+iadg];
    jh = [jh1;jadg];
    sh = [sh1;sadg];
    h = sparse(ih,jh,sh, 2*n*m, 1);
    %clearvars -except h H
    z = H\h;
    y = z(1:n*m);
    p = z(n*m+1:2*n*m);

    
    u(:,1) = Einaktiv1*(-V1*p-v)+uaktiv(:,1); 
    u(:,2) = Einaktiv2*(-V2*p-v)+uaktiv(:,2); 
    u(:,3) = Einaktiv3*(-V3*p-v)+uaktiv(:,3); 
    u(:,4) = Einaktiv4*(-V4*p-v)+uaktiv(:,4); 
 
yoptvec = y;
poptvec = p;
for i = 1:numcontrols
    uopt(:,i) = min(ub(:,i),max(ua(:,i),u(:,i))); 
end