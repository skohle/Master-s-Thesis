function  [ypod] = stateeqn_pod(deltat, y0pod, c, fxyt, au, ui, Mpod , Kpod , POD , M , Mpod4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Sophia Kohle
% -------------------------------------------------------------------------
%
% Computing the state y to a given control u
%
% y_t(x,t) - /\y(x,t) + c y^3(x,t) + d(x,t) = sum(i=1:4)au_i(x) u_i(t) 
% y(x,0) = y_0(x) 
% dy/dnu(x,t) = 0
%
% -------------------------------------------------------------------------
% 
% Used algorithms:  POD discretization in space, semi-implicit 
%                   Euler method in time
% 
% -------------------------------------------------------------------------
%
% INPUT:
%
% deltat    % time mesh size
% y0pod     % initial condition POD
% c
% fxyt      % parameter f(x,t) in PDE 
% au        % parameter au in PDE   
% ui        % given control
% Mpod,Kpod % POD matrices
% POD       % POD basis
% M         % massmatrix FEM
% Mpod4     % POD matrix phi^4
% 
% -------------------------------------------------------------------------
%
% OUTPUT:
%
% ypod              % POD solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[n,numcontrols] = size(ui);
[m,r] = size(POD);


% solve the PDE
% ------------------------------------------------------------------------

[iy,jy,sy] = find(y0pod);
ypre = sparse(iy,jy,sy,r,1);

% compute RHS
for k = 2:n    
    Muk = 0;    
    for i = 1:numcontrols
        Muk = Muk + POD'*M*au(:,i)*ui(k,i);
    end 
%     
    
% Handle the nonlinear term with POD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:r
    for j=1:r
        My(i,j) = ypre'*Mpod4(:,:,i,j)*ypre;
    end
end


%semi-implicit Euler
  [i,j,s] = find((Mpod/deltat + Kpod)\((Mpod*ypre)/deltat + Muk - c*My*ypre- POD'*M*fxyt(:,k)));
    iy = [iy,i];
    jy = [jy,j+k-1];
    sy = [sy;s];
    ypre = sparse(i,j,s,r,1);
end
ypod = sparse(iy,jy,sy,r,n);

clear iy jy sy i s;






