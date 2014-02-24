function [yopt,uopt,popt] = ssn_fem(n,m,xy,tri,deltat,lambda,c,dxyt,ay,au,ua,ub,ustart,ystart,pstart,M,K,y0,yq,V1,V2,V3,V4,v,Mt)

% =========================================================================
% Author: Sophia Kohle, Technische Universität Berlin
% =========================================================================
%
% Task: Computation of an optimal control for the problem
%
% min 1/2 ||y-yq||^2 + lambda/2 sum(i=1:4) ||u_i||^2 + (ay,y(T))
% s.t.
%       y_t(x,t) - /\y(x,t) + c y^3(x,t) + d(x,t) = sum(i=1:4)au_i(x) u_i(t) 
%       y(x,0) = y_0(x) 
%       dy/dnu(x,t) = 0
%
%       -1 <= u_i(t) <= 1     , i = 1:4
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Numerical methods:    Semi smooth Newton Method, Finite Element Method,
%                       Primal-dual Active Set Method
%              
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Input parameters:
%
% m             % number of finite elements in [0,L]
% n             % number of time instances in [0,T]
% xy, tri       % points, triangles of FEM dicretization
% lambda        % Tikhonov parameter (small)
% yq            % term yq of objective functional
% ay            % term ay of objective functional
% au            % term au of objective functional
% dxyt          % term dxyt on left-hand side of PDE 
% y0            % initial value
% ua            % lower control bound
% ub            % upper control bound
% c             % parameter c in PDE
% ystart        % start value y^0 in SQP method
% ustart        % start value u^0 in SQP method
% pstart        % start value p^0 in SQP method
% M,K           % FE matrices
% V1,V2,V3,V4,v % matrices for projection term
% deltat        % step size in time
% numcontrols   % number of control functions
% Mt        
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Output:   optimal solution set (uopt, yopt,popt)
% =========================================================================

% Parameters
TOL = 1e-3;
countermax = 30;


yk = ystart;
uk = ustart;
pk = pstart;
diff = Inf;
counter = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Semi smooth Newton
%   psi(y)
%   ==>
%   psi(yk) + psi'(yk)(y-yk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (diff > TOL && counter < countermax)
    
    ak = 3*c*(yk.^2);
    bk = 2*c*(yk.^3) - dxyt;
    ck = ak;
    dk = - 6*c*(yk.*pk);
    ek = 6*c*((yk.^2).*pk) - yq;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Discretization of state equation in time and space
    % G1 y  = sum(G2i ui) + g
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [iI,jI,sI] = find(speye(m));
    [iM,jM,sM] = find(M);
    [iy0,jy0,sy0] = find(y0);

    % G1(1:m,1:m) = I;
    iG1 = iI;
    jG1 = jI;
    sG1 = sI;
    % g(1:m) = y0;
    ig = iy0;
    jg = jy0;
    sg = sy0;

    for k = 1:(n-1)
        % G1((k*m+1):((k+1)*m),((k-1)*m+1):(k*m)) = -M/deltat;
        iG1 = [iG1;iM+k*m];
        jG1 = [jG1;jM+(k-1)*m];
        sG1 = [sG1;-sM./deltat];
        % G1((k*m+1):((k+1)*m),(k*m+1):((k+1)*m)) = M/deltat + K + diag(M*akvec(:,k+1));
%             aktri(:,1) = pdeintrp(xy,tri,ak(:,k+1));
%             [~,Mak,~] = assema(xy,tri,0,aktri',0);
%             [iMKak,jMKak,sMKak] = find(M/deltat + K + Mak);
        % empirical better?
        [iMKak,jMKak,sMKak] = find(M/deltat + K + diag(M*ak(:,k+1)));
        iG1 = [iG1;iMKak+k*m];
        jG1 = [jG1;jMKak+k*m];
        sG1 = [sG1;sMKak];
        % g((k*m+1):((k+1)*m),1)                 = M*bkvec(:,k+1)
        [iMbk,jMbk,sMbk] = find(M*bk(:,k+1));
        ig = [ig;iMbk+k*m];
        jg = [jg;jMbk];
        sg = [sg;sMbk];
    end
    G1 = sparse(iG1,jG1,sG1,n*m,n*m);
    g = sparse(ig,jg,sg,n*m,1);

    % G21
    iG2 = [];
    jG2 = [];
    sG2 = [];
    [iMaui,jMaui,sMaui] = find(M*au(:,1));
    for k = 1:(n-1)
        % G21((k*m+1):((k+1)*m),k+1)             = M*au(:,1);
        iG2 = [iG2;iMaui+k*m];
        jG2 = [jG2;jMaui+k];
        sG2 = [sG2;sMaui];
    end
    G21 = sparse(iG2,jG2,sG2,n*m,n);
    % G22
    iG2 = [];
    jG2 = [];
    sG2 = [];
    [iMaui,jMaui,sMaui] = find(M*au(:,2));
    for k = 1:(n-1)
        % G22((k*m+1):((k+1)*m),k+1)             = M*au(:,2);
        iG2 = [iG2;iMaui+k*m];
        jG2 = [jG2;jMaui+k];
        sG2 = [sG2;sMaui];
    end
    G22 = sparse(iG2,jG2,sG2,n*m,n);        
    % G23
    iG2 = [];
    jG2 = [];
    sG2 = [];
    [iMaui,jMaui,sMaui] = find(M*au(:,3));
    for k = 1:(n-1)
        % G23((k*m+1):((k+1)*m),k+1)             = M*au(:,3);
        iG2 = [iG2;iMaui+k*m];
        jG2 = [jG2;jMaui+k];
        sG2 = [sG2;sMaui];
    end
    G23 = sparse(iG2,jG2,sG2,n*m,n);        
    % G24
    iG2 = [];
    jG2 = [];
    sG2 = [];
    [iMaui,jMaui,sMaui] = find(M*au(:,4));
    for k = 1:(n-1)
        % G24((k*m+1):((k+1)*m),k+1)             = M*au(:,4);
        iG2 = [iG2;iMaui+k*m];
        jG2 = [jG2;jMaui+k];
        sG2 = [sG2;sMaui];
    end
    G24 = sparse(iG2,jG2,sG2,n*m,n);        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Discretization of adjoint equation in time and space
    % ADG1 p = ADG2 y + adg
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [iay,jay,say] = find(ay);
    % ADG1((n-1)*m+1:n*m,(n-1)*m+1:n*m) = I;
    iADG1 = iI+(n-1)*m;
    jADG1 = jI+(n-1)*m;
    sADG1 = sI;
    % ADG2
    iADG2 = [];
    jADG2 = [];
    sADG2 = [];
    % adg((n-1)*m+1:n*m) = ay;
    iadg = iay+(n-1)*m;
    jadg = jay;
    sadg = say;

    for k = 1:(n-1)
        % ADG1(((k-1)*m+1):(k*m),((k-1)*m+1):(k*m)) = M/deltat + K + Mck(k);
        cktri(:,1) = pdeintrp(xy,tri,ck(:,k));
        [~,Mck,~] = assema(xy,tri,0,cktri',0);
        [iMKck,jMKck,sMKck] = find(M/deltat + K + Mck);
        iADG1 = [iADG1;iMKck+(k-1)*m];
        jADG1 = [jADG1;jMKck+(k-1)*m];
        sADG1 = [sADG1;sMKck];
        % ADG1(((k-1)*m+1):(k*m),(k*m+1):((k+1)*m)) = - M/deltat;
        iADG1 = [iADG1;iM+(k-1)*m];
        jADG1 = [jADG1;jM+k*m];
        sADG1 = [sADG1;-sM./deltat];
        % ADG2(((k-1)*m+1):(k*m),((k-1)*m+1):(k*m)) = M+Mdk(t_k);
        dktri(:,1) = pdeintrp(xy,tri,dk(:,k));
        [~,Mdk,~] = assema(xy,tri,0,dktri',0);
        [iD,jD,sD] = find(M+Mdk);
        iADG2 = [iADG2;iD+(k-1)*m];
        jADG2 = [jADG2;jD+(k-1)*m];
        sADG2 = [sADG2;sD];
        % adg(((k-1)*m+1):(k*m)) = M*ek(t_k);
        [iek,jek,sek] = find(M*ek(:,k));
        iadg = [iadg;iek+(k-1)*m];
        jadg = [jadg;jek];
        sadg = [sadg;sek];
    end

    ADG1 = sparse(iADG1,jADG1,sADG1,n*m,n*m);
    ADG2 = sparse(iADG2,jADG2,sADG2,n*m,n*m);
     adg = sparse(iadg,jadg,sadg,n*m,1);        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Primal Dual Active Set Method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [yl,ul,pl] = active_set(ua,ub,G1,G21,G22,G23,G24,g,ADG1,ADG2,adg,V1,V2,V3,V4,v,m,n,lambda,reshape(yk,n*m,1),reshape(pk,n*m,1));

    yl = reshape(yl,m,n);
    pl = reshape(pl,m,n);


    % distance between last two solutions
    diffy = error_numerical_Q(yl,yk,M,Mt);
    for i = 1:4
        diffu(i) = (ul(:,i)-uk(:,i))'*Mt*(ul(:,i)-uk(:,i));
    end
    diff = diffy + sum(diffu);

    yk = yl;
    uk = ul;
    pk = pl;

    counter = counter+1;


end % SSN
disp(' ');
disp(['Number of iterations: ' num2str(counter)]);
disp(['Stopping criteria  ||ul-uk||^2+||yl-yk||^2 = ' num2str(diff)]);

yopt = yl;
uopt = ul;
popt = pl;



