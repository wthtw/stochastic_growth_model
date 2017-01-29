%% Solving a stochastic growth model with collocation method
tic;

clc;clear;

sigma = 1.50;       % utility parameter
delta = 0.10;       % depreciation rate
beta = 0.95;        % discount factor
alpha = 0.30;       % capital elasticity of output
rho = 0.80;         % persistence of the shock on Z
epsi = 1e-5;       % convergence parameter

% Discretization of the shocks, Markov AR(1) process
nbe = 15;           % number of shocks, 15
se = 0.12;          % volatility of the shock
mu = 0;             % mean of shocks residuals iid normal process

% Matrice of transition probabilities using tauchen and hussey
muZ = 0.5;            % mean of shocks
[Z,P] = tauchenhussey(nbe,muZ,rho,se,se);

% constructing grid values for K, using collocation method
nbk     = 50;        % number of collocation points, 50
kmin  = 0.2;
kmax  = 6;
basis = fundefn('lin',nbk,kmin,kmax);  % cubic spline basis   
kgrid = funnode(basis);                    
Phi   = funbas(basis);

% setting initial guess for basis coefficients c
c = zeros(nbk,nbe);
%load c;
%load c_gscoll;

% initialization of various matrices
Tv = zeros(nbk,nbe); 
dr = zeros(nbk,nbe);
C = zeros(nbk,nbe);
util = zeros(nbk,nbe);

% starting iterations
for it=1:1000
    
   cold   = c;
      
   for i=1:nbk   
        for k = 1:nbe
        C(:,k) = Z(k)*(kgrid(i)^alpha) + (1-delta)*kgrid(i) - kgrid;
        util(:,k) = (C(:,k).^(1-sigma)-1)/(1-sigma);
        end

        neg = C<0;
        util(neg) = -1e12;
        [Tv(i,:), dr(i,:)] = max(util + beta*(funeval(c,basis,kgrid)*P));

   end

       c = Phi\Tv;
       crit = norm(c-cold);
       if crit<epsi, break,end
              
end

Kp = kgrid(dr); % getting the optimal policy K'* for each K in the grid

% computing the optimal policy C* for each K in the grid
for k=1:nbe;
C(:,k) = Z(k).*(kgrid.^alpha) + (1-delta).*kgrid - Kp(:,k);
end
neg = C<0;
C(neg) = 1e-6;
util(neg) = -1e12;

t = toc; % getting the duration in processor time

save ('Kp','Kp');
save ('C','C');
save ('t','t');
