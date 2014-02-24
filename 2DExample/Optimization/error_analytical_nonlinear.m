function [error,mineig,cpu_hess] = error_analytical_nonlinear(ui,y,p,ua,ub,c,dxyt,y0,au,M,K,lambda,deltat,V1,V2,V3,V4,v1,Mt,xy,tri)

% =========================================================================
% Author: Eileen Kammann, Technische Universität Berlin
% =========================================================================
%
%
% Task: A-posteriori error estimation for parabolic nonlinear problems like
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
%   ||uopt - uoptpod||_L²(0,T;R4) <= ||sigma||_L²(0,T;R4)/mineig*deltat
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Input parameters:
%
% ui            % controls ui
% y             % state y
% p             % adjoint state p
% ua            % lower control bound
% ub            % upper control bound
% xy, tri       % points, triangles of FEM dicretization
% lambda        % Tikhonov parameter
% au            % term au of objective functional
% dxyt          % term dxyt on left-hand side of PDE 
% y0            % initial value
% c             % parameter c in PDE
% M,K           % FE matrices
% V1,V2,V3,V4,v % matrices for projection term
% deltat        % step size in time
% Mt            % time discr. matrix
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Output:   analytical error estimate
% =========================================================================


[m,n] = size(p);
numcontrols = size(ui,2);
p = reshape(p,m*n,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) perturbation vector sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vi = (lambda ui + (aui,p)_L²(Omega) )
v = zeros(n,numcontrols);
v(:,1) = lambda*(ui(:,1) + V1*p + v1);
v(:,2) = lambda*(ui(:,2) + V2*p + v1);
v(:,3) = lambda*(ui(:,3) + V3*p + v1);
v(:,4) = lambda*(ui(:,4) + V4*p + v1);

% sigmai  = (lambda ui + (aui,p)_L²(Omega) )- , if u = u_a
%           (lambda ui + (aui,p)_L²(Omega) )+ , if u = u_b
%          -(lambda ui + (aui,p)_L²(Omega) )  , else
sigma = zeros(n,numcontrols);
for i = 1:numcontrols
    for k = 1:n
        if ui(k,i) == ua(k,i)
            if v(k,i) >= 0
                sigma(k,i) = 0;
            else
                sigma(k,i) = -v(k,i);
            end
        elseif ui(k,i) == ub(k,i)
            if v(k,i) <= 0
                sigma(k,i) = 0;
            else
                sigma(k,i) = v(k,i);
            end
        else
            sigma(k,i) = -v(k,i);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) smallest eigenvalue of reduced hessian matrix 
%     (approximated by last SQP Hessian)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = reshape(p,m,n);

        timestart = cputime;

    ak = 3*c*(y.^2);
    bk = 2*c*(y.^3) - dxyt;
    dk = -6*c*(y.*p);


    % state equation
    % G1 y  = G2 u + g
    % --------------------------------------------------------------------
    % G1(1:m,1:m) = M;
    % g(1:m) = M*y0;
    % for k = 1:(n-1)
    %      G1((k*m+1):((k+1)*m),((k-1)*m+1):(k*m)) = -M/deltat;
    %      G1((k*m+1):((k+1)*m),(k*m+1):((k+1)*m)) = M/deltat + K + diag(M*ak(k+1));
    %      G2i((k*m+1):((k+1)*m),k+1)              = Maui;
    %       g((k*m+1):((k+1)*m),1)                 = M*bk(k+1)
    % end
    [iM,jM,sM] = find(M);
    [iy0,jy0,sy0] = find(M*y0);

    % G1(1:m,1:m) = M;
    iG1 = iM;
    jG1 = jM;
    sG1 = sM;
    % g(1:m) = y0vec;
    ig = iy0;
    jg = jy0;
    sg = sy0;

    for k = 1:(n-1)
        % G1((k*m+1):((k+1)*m),((k-1)*m+1):(k*m)) = M ./ (-deltat);
        iG1 = [iG1;iM+k*m];
        jG1 = [jG1;jM+(k-1)*m];
        sG1 = [sG1;-sM./deltat];
        % G1((k*m+1):((k+1)*m),(k*m+1):((k+1)*m)) = M ./ deltat + K + diag(M*akvec(:,k+1));
%             aktri(:,1) = pdeintrp(xy,tri,akvec(:,k+1));
%             [~,Mak,~] = assema(xy,tri,0,aktri',0);
%             [iMKak,jMKak,sMKak] = find(M./deltat + K + Mak);
        % empirical better
        [iMKak,jMKak,sMKak] = find(M./deltat + K + diag(M*ak(:,k+1)));
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
    
    
    G2 = [G21,G22,G23,G24];

    % null space matrix
    Z(1:n*m,1:n*numcontrols) = G1 \ G2;
    Z(n*m+1:n*m+n*numcontrols,1:n*numcontrols) = speye(numcontrols*n);

    [iMt,jMt,sMt] = find(Mt);

    iQ = [];
    jQ = [];
    sQ = [];
    for k = 1:n
        % ||y² - 6c*yp*pp*y²||_L²(Q)
                dktri(:,1) = pdeintrp(xy,tri,dk(:,k));
                [~,Mdk,~] = assema(xy,tri,0,dktri',0);
                [iMdk,jMdk,sMdk] = find(deltat*(M + Mdk));
                if k == 1 || k == n
                    iQ = [iQ;iMdk+(k-1)*m];
                    jQ = [jQ;jMdk+(k-1)*m];
                    sQ = [sQ;sMdk/2];
                else
                    iQ = [iQ;iMdk+(k-1)*m];
                    jQ = [jQ;jMdk+(k-1)*m];
                    sQ = [sQ;sMdk];
                end
    end
    for i = 1:numcontrols
        % ||lambda*ui||_L²(0,T)
        iQ = [iQ;iMt+m*n+(i-1)*n];
        jQ = [jQ;jMt+m*n+(i-1)*n];
        sQ = [sQ;lambda*sMt];
    end

    Q = sparse(iQ,jQ,sQ,n*m+n*numcontrols,n*m+n*numcontrols);

    % reduced Hessian
    Hess = Z'*Q*Z;
    Hess = (Hess + Hess')/2;
        cpu_hess = cputime - timestart;    
    eigHess = eigs(Hess);
    mineig = min(eigHess);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3)  ||uopt - uoptpod||_L²(0,T;R4) <= ||sigma||_L²(0,T;R4)/mineig*deltat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normsigmasquared = zeros(numcontrols,1); 
for i = 1:numcontrols
    normsigmasquared(i) = sigma(:,i)'*Mt*sigma(:,i);
end
error = sqrt(sum(normsigmasquared))/mineig*deltat;