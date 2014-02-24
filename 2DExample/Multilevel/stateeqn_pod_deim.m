function  ypoddeim = stateeqn_pod_deim(deltat, y0pod, c, fxyt, au, uiopt, Mpod , Kpod , POD ,U,P, M )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Sophia Kohle
% -------------------------------------------------------------------------
%
% Computing the state y to a given control uiopt
%
% y_t(x,t) - /\y(x,t) + c y^3(x,t) + d(x,t) = sum(i=1:4)au_i(x) u_i(t) 
% y(x,0) = y_0(x) 
% dy/dnu(x,t) = 0
%
% -------------------------------------------------------------------------
% 
% Used algorithms:  POD-DEIM discretization in space, semi-implicit 
%                   Euler method in time
% 
% -------------------------------------------------------------------------
%
% INPUT:
%
% deltat        % time mesh size
% y0pod         % initial POD condition
% c
% fxyt          % parameter f(x,t) in PDE 
% au            % parameter au in PDE   
% uiopt         % given POD control 
% Mpod , Kpod   % POD matrices
% POD           % POD basis
% U,P           % POD-DEIM matrices
% M             % massmatrix FEM
% x1x2          % discretization of omega
% r_deim        % number of DEIM-interpolation points
% -------------------------------------------------------------------------
%
% OUTPUT:
%
% ypoddeim  - optimal state computed with POD-DEIM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,numcontrols] = size(uiopt);
[m,r] = size(POD);

% solve the PDE
% ------------------------------------------------------------------------

[iy,jy,sy] = find(y0pod);
ypre = sparse(iy,jy,sy,r,1);

% Handle the nonlinear term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Precomputations
Proj=(POD'*M*U)/(P'*U);

PPOD=P'*POD;

PODM=POD'*M;
% compute RHS
for k = 2:n
    Muk = 0;  
    for i = 1:numcontrols
        Muk = Muk + PODM*au(:,i)*uiopt(k,i);
    end 

%semi-implicit Euler with POD-DEIM
  [i,j,s] = find((Mpod/deltat + Kpod)\((Mpod*ypre)/deltat + Muk -c*Proj*(PPOD*ypre).^3- PODM*fxyt(:,k)));
    iy = [iy,i];
    jy = [jy,j+k-1];
    sy = [sy;s];
    ypre = sparse(i,j,s,r,1);
end
ypoddeim = sparse(iy,jy,sy,r,n);

clear iy jy sy i s;






