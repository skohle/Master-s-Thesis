function [yopt,uopt,popt] = ssn_pod_deim(n,r,xy,tri,deltat,lambda,c,dxyt,ay,au,ua,ub,ustart,ystart,pstart,M,K,Mt,y0,yq,POD,U,P,Mfem)

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
%       ua <= u_i(t) <= ub     , i = 1:4
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Numerical methods:    SSN Method, POD Method,
%                       Primal-dual Active Set Method
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Input parameters:
%
% n         % number of time instances in [0,T]
% r         % rank of POD basis
% lambda    % Tikhonov parameter (small)
% yd        % term yd of objective functional
% ay        % term ay of objective functional
% au        % term au of objective functional
% bt        % term bt on right-hand side of PDE 
% y0        % initial value
% ua        % lower control bound
% ub        % upper control bound
% delta     % parameter delta in PDE   
% c         % parameter c in PDE
% usnap     % control for generation of snapshots
% ystart    % start value y^0 in SQP method
% ustart    % start value u^0 in SQP method
% deltax    % step size in space            
% M,K,Q,G   % FE matrices
% deltat    % step size in time
% Mt        
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                
% Output:   optimal solution set (uopt, yopt,popt)
% =========================================================================

% Parameters
countermax = 30;
TOL = 1e-7;
numcontrols = size(au,2);
yk = ystart;
uk = ustart;
pk = pstart;
diff = Inf;
diffu = zeros(numcontrols,1);
counter = 0;

% interpolate POD to triangles
for i = 1:r
    PODtri(:,i) = pdeintrp(xy,tri,POD(:,i));
end

numtri = size(PODtri,1);
Atri = zeros(numtri,1);
for i = 1:numtri
    x1 = xy(1,tri(1,i));
    y1 = xy(2,tri(1,i));
    x2 = xy(1,tri(2,i));
    y2 = xy(2,tri(2,i));
    x3 = xy(1,tri(3,i));
    y3 = xy(2,tri(3,i));
    Atri(i) = 0.5*((x2-x1)*(y3-y1) - (x3-x2)*(y2-y1));
end
    
% integral over Omega of product of each four POD basis functions
Mpod4 = zeros(r,r);
for ip = 1:r
    for iq = 1:r
        for i = 1:r
            for j = 1:r
                Mpod4(i,j,ip,iq) = sum(Atri.*PODtri(:,i).*PODtri(:,j).*PODtri(:,ip).*PODtri(:,iq));
            end
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection term from variational inequality
% ui(t) = -1/lambda int_Omega( au_i(x)p(x,t)dx )
% ==> ui = - Vi*p - v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V1
[iau,jau,sau] = find((au(:,1)'*M)./lambda);
iV = [];
jV = [];
sV = [];
for k = 1:n
    iV = [iV,k*iau];
    jV = [jV,(k-1)*r+jau];
    sV = [sV,sau];
end
V1 = sparse(iV,jV,sV,n,n*r);
% V2
[iau,jau,sau] = find((M*au(:,2))'./lambda);
iV = [];
jV = [];
sV = [];
for k = 1:n
    iV = [iV,k*iau];
    jV = [jV,(k-1)*r+jau];
    sV = [sV,sau];
end
V2 = sparse(iV,jV,sV,n,n*r);
% V3
[iau,jau,sau] = find((M*au(:,3))'./lambda);
iV = [];
jV = [];
sV = [];
for k = 1:n
    iV = [iV,k*iau];
    jV = [jV,(k-1)*r+jau];
    sV = [sV,sau];
end
V3 = sparse(iV,jV,sV,n,n*r);
% V4
[iau,jau,sau] = find((M*au(:,4))'./lambda);
iV = [];
jV = [];
sV = [];
for k = 1:n
    iV = [iV,k*iau];
    jV = [jV,(k-1)*r+jau];
    sV = [sV,sau];
end
V4 = sparse(iV,jV,sV,n,n*r);

v = sparse(n,1);

Proj=sparse((POD'*Mfem*U)/(P'*U));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SSN
%   psi(y)
%   ==>
%   psi(yk) + psi'(yk)(y-yk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (diff > TOL && counter < countermax)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Discretization of state equation in time and space
    % G1 y  = G2 u + g
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Discretization of adjoint equation in time and space
    % ADG1 p = ADG2 y + adg
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [iM,jM,sM] = find(M);
    [iI,jI,sI] = find(speye(r));
    [iy0,jy0,sy0] = find(y0);
    [iay,jay,say] = find(ay);

    % G1(1:m,1:m) = I;
    iG1 = iI;
    jG1 = jI;
    sG1 = sI;
    % g(1:m) = y0;
    ig = iy0;
    jg = jy0;
    sg = sy0;

    % ADG1((n-1)*m+1:n*m,(n-1)*m+1:n*m) = I;
    iADG1 = iI+(n-1)*r;
    jADG1 = jI+(n-1)*r;
    sADG1 = sI;
    % ADG2
    iADG2 = [];
    jADG2 = [];
    sADG2 = [];
    % adg((n-1)*m+1:n*m) = ay;
    iadg = iay+(n-1)*r;
    jadg = jay;
    sadg = say;
        
    for k = 1:(n-1)
        
        % ADG1(((k-1)*m+1):(k*m),((k-1)*m+1):(k*m)) = M/deltat + K + Mck(k);
        % ck = 3*c*(yk.^2);
        if k == 1
            Mck = diag(Proj*P'*(POD*yk(:,k)).^2);
        else
            Mck = Mak;
        end
        [iMKck,jMKck,sMKck] = find(M/deltat + K + 3*c*Mck);
        iADG1 = [iADG1;iMKck+(k-1)*r];
        jADG1 = [jADG1;jMKck+(k-1)*r];
        sADG1 = [sADG1;sMKck];
        % ADG1(((k-1)*m+1):(k*m),(k*m+1):((k+1)*m)) = - M/deltat;
        iADG1 = [iADG1;iM+(k-1)*r];
        jADG1 = [jADG1;jM+k*r];
        sADG1 = [sADG1;-sM/deltat];
        % ADG2(((k-1)*m+1):(k*m),((k-1)*m+1):(k*m)) = M+Mdk(t_k);
        % dk = - 6*c*(yk.*pk);
        Mdk = zeros(r,r);
        for ip = 1:r
            for iq = 1:r
                Mdk(ip,iq) = -6*c*yk(:,k)'*Mpod4(:,:,ip,iq)*pk(:,k);
            end
        end
        [iD,jD,sD] = find(M + Mdk);
        iADG2 = [iADG2;iD+(k-1)*r];
        jADG2 = [jADG2;jD+(k-1)*r];
        sADG2 = [sADG2;sD];
        % adg(((k-1)*m+1):(k*m))                    = M*ek(t_k);
        % ek = 6*c*((yk.^2).*pk) - yq;
        [iek,jek,sek] = find(6*c*Mck*pk(:,k) - M*yq(:,k));
        iadg = [iadg;iek+(k-1)*r];
        jadg = [jadg;jek];
        sadg = [sadg;sek];
        
        % G1((k*m+1):((k+1)*m),((k-1)*m+1):(k*m)) = -M/deltat;
        iG1 = [iG1;iM+k*r];
        jG1 = [jG1;jM+(k-1)*r];
        sG1 = [sG1;-sM/deltat];
        % ak = 3*c*(yk.^2);
        Mak = diag(Proj*P'*(POD*yk(:,k+1)).^2);
        % G1((k*m+1):((k+1)*m),(k*m+1):((k+1)*m)) = M/deltat + K + Mak;
        [iMKak,jMKak,sMKak] = find(M/deltat + K + 3*c*Mak);
        iG1 = [iG1;iMKak+k*r];
        jG1 = [jG1;jMKak+k*r];
        sG1 = [sG1;sMKak];
        % g((k*m+1):((k+1)*m),1) = 2c Mak*yk(:,k+1) - M*dyxt(:,k+1)
        % bk = 2*c*(yk.^3) - dxyt;
        [iMbk,jMbk,sMbk] = find(2*c*Mak*yk(:,k+1) - M*dxyt(:,k+1));
        ig = [ig;iMbk+k*r];
        jg = [jg;jMbk];
        sg = [sg;sMbk];
        
    end
    G1 = sparse(iG1,jG1,sG1,n*r,n*r);
    g = sparse(ig,jg,sg,n*r,1);

    % G21
    iG2 = [];
    jG2 = [];
    sG2 = [];
    [iMaui,jMaui,sMaui] = find(M*au(:,1));
    for k = 1:(n-1)
        % G21((k*m+1):((k+1)*m),k+1) = M*au(:,1);
        iG2 = [iG2;iMaui+k*r];
        jG2 = [jG2;jMaui+k];
        sG2 = [sG2;sMaui];
    end
    G21 = sparse(iG2,jG2,sG2,n*r,n);
    % G22
    iG2 = [];
    jG2 = [];
    sG2 = [];
    [iMaui,jMaui,sMaui] = find(M*au(:,2));
    for k = 1:(n-1)
        % G22((k*m+1):((k+1)*m),k+1) = M*au(:,2);
        iG2 = [iG2;iMaui+k*r];
        jG2 = [jG2;jMaui+k];
        sG2 = [sG2;sMaui];
    end
    G22 = sparse(iG2,jG2,sG2,n*r,n);        
    % G23
    iG2 = [];
    jG2 = [];
    sG2 = [];
    [iMaui,jMaui,sMaui] = find(M*au(:,3));
    for k = 1:(n-1)
        % G23((k*m+1):((k+1)*m),k+1) = M*au(:,3);
        iG2 = [iG2;iMaui+k*r];
        jG2 = [jG2;jMaui+k];
        sG2 = [sG2;sMaui];
    end
    G23 = sparse(iG2,jG2,sG2,n*r,n);        
    % G24
    iG2 = [];
    jG2 = [];
    sG2 = [];
    [iMaui,jMaui,sMaui] = find(M*au(:,4));
    for k = 1:(n-1)
        % G24((k*m+1):((k+1)*m),k+1) = M*au(:,4);
        iG2 = [iG2;iMaui+k*r];
        jG2 = [jG2;jMaui+k];
        sG2 = [sG2;sMaui];
    end
    G24 = sparse(iG2,jG2,sG2,n*r,n);        

    ADG1 = sparse(iADG1,jADG1,sADG1,n*r,n*r);
    ADG2 = sparse(iADG2,jADG2,sADG2,n*r,n*r);
     adg = sparse(iadg,jadg,sadg,n*r,1);        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Primal Dual Active Set Method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [yl,ul,pl] = active_set(ua,ub,G1,G21,G22,G23,G24,g,ADG1,ADG2,adg,V1,V2,V3,V4,v,r,n,lambda,reshape(yk,n*r,1),reshape(pk,n*r,1));
    
    yl = reshape(yl,r,n);
    pl = reshape(pl,r,n);


    % distance between last two solutions
    diffy = error_numerical_Q(yl,yk,M,deltat);
    for i = 1:numcontrols
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
% Backtransformation
yopt = POD*yl;
uopt = ul;
popt = POD*pl;

