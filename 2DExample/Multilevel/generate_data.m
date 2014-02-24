function generate_data( mmin , n )

% =========================================================================
% Author: Sophia Kohle, Technische Universit�t Berlin
% =========================================================================
%
% Masterthesis
% 2D Example_LS
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Task: Computation of an optimal control for the problem
%
% min 1/2 ||y-yq||^2 + lambda/2 sum(i=1:4) ||u_i||^2 + (ay,y(T))
% s.t.
%       y_t(x,t) - /\y(x,t) + c * y^3(x,t) + f(x,t) = sum(i=1:4)au_i(x) u_i(t) 
%       y(x,0) = y_0(x) 
%       dy/dnu(x,t) = 0
%
%       ua <= u_i(t) <= ub     , i = 1:4
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% 
% Generation of the mesh and the needed functions
% 


disp(' ');
disp('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*');

% Parameters (1)
% ------------------------------------------------------------------------
omegashape = rectangle_g(0,pi,0,pi);          % geometry file for space Omega
gammashape = rectangle_b('0', '0', '0', '0'); % boundary file for space Omega
T = 1;                       % time interval (0,T)
lambda = 0.01;               % Tikhonov parameter (small)
c = 1;                       % parameter c in PDE
numcontrols = 4;             % number of control functions
ua = -ones(n,numcontrols);   % lower control bound
ub = ones(n,numcontrols);    % upper control bound

%TOL1 = 1*1e-5;

s0=5;
reset=5;

% Discretization of (0,T) piecewise constant
deltat = T/(n-1);
t = (0:deltat:T)';
Mt = deltat * eye(n,n);
Mt(1,1) = Mt(1,1)/2;
Mt(n,n) = Mt(n,n)/2;
   
% Discretization of Omega by Finite Element method
[m,x1x2,edges,tri] = grid_generation(mmin, omegashape);
disp(' ');
disp(['   Number of Finite Elements: m = ' num2str(m)]); 
disp(' ');
x1 = x1x2(1,:)';

x2 = x1x2(2,:)';
au(:,1) = max(0,10-50*(x1-pi/4).^2-50*(x2-pi/4).^2);      % term au of objective functional
au(:,2) = max(0,10-50*(x1-3/4*pi).^2-50*(x2-pi/4).^2);
au(:,3) = max(0,10-50*(x1-3/4*pi).^2-50*(x2-3/4*pi).^2);
au(:,4) = max(0,10-50*(x1-pi/4).^2-50*(x2-3/4*pi).^2);
[M,K,~,~,~,~,~] = assem_FEM(x1x2, edges, tri, gammashape, au);

% Parameters (2)
% ------------------------------------------------------------------------
ay = cos(x1).*cos(x2)*lambda*T^2;	% term ay of objective functional
y0 = cos(x1).*cos(x2);            % initial value      

yq = zeros(m,n);                % term yq of objective functional
for i = 1:m
    for j = 1:n
        yq(i,j) = cos(x1(i))*cos(x2(i))*(1 + 2*lambda*t(j) - 2*lambda*t(j)^2 - 3*c*lambda*t(j)^2*cos(x1(i))^2*cos(x2(i))^2);
    end
end


% constructed optimal solutions
% ------------------------------------------------------------------------
 
    % optimal control
        pxy = cos(x1).*cos(x2);
        Aaup = au'*M*pxy;
    uiopt = zeros(n,numcontrols);
    for i = 1:numcontrols
        uiopt(:,i) = min(ub(:,i),max(ua(:,i),-t.^2*Aaup(i)));
    end

        uopt = zeros(m,n);
        for i = 1:n
            for j = 1:m
                for k = 1:numcontrols
                    uopt(j,i) = uopt(j,i) + au(j,k)*uiopt(i,k);
                end
            end
        end

% Parameters (3)
fxyt = zeros(m,n);              % term f(x,y,t) on left-hand side of PDE 
for i = 1:n
    fxyt(:,i) = uopt(:,i) - 2*cos(x1).*cos(x2) - c*cos(x1).^3.*cos(x2).^3;
end
fxyt = sparse(fxyt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection term from variational inequality
% ui(t) = -1/lambda int_Omega( au_i(x)p(x,t)dx )
% ==> ui = - Vi*p - v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V1
[iF,jF,sF] = find((au(:,1)'*M)./lambda);
iV = [];
jV = [];
sV = [];
for k = 1:n
    iV = [iV,k*iF];
    jV = [jV,(k-1)*m+jF];
    sV = [sV,sF];
end
V1 = sparse(iV,jV,sV,n,n*m);
% V2
[iF,jF,sF] = find((M*au(:,2))'./lambda);
iV = [];
jV = [];
sV = [];
for k = 1:n
    iV = [iV,k*iF];
    jV = [jV,(k-1)*m+jF];
    sV = [sV,sF];
end
V2 = sparse(iV,jV,sV,n,n*m);
% V3
[iF,jF,sF] = find((M*au(:,3))'./lambda);
iV = [];
jV = [];
sV = [];
for k = 1:n
    iV = [iV,k*iF];
    jV = [jV,(k-1)*m+jF];
    sV = [sV,sF];
end
V3 = sparse(iV,jV,sV,n,n*m);
% V4
[iF,jF,sF] = find((M*au(:,4))'./lambda);
iV = [];
jV = [];
sV = [];
for k = 1:n
    iV = [iV,k*iF];
    jV = [jV,(k-1)*m+jF];
    sV = [sV,sF];
end
V4 = sparse(iV,jV,sV,n,n*m);

v = sparse(n,1);

save data.mat