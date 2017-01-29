%% Solving a stochastic growth model with Function iteration
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

%  constructing grid values for K, evenly spaced
nbk     = 50; % number of points, 50
kmin  = 0.2;   % min
kmax  = 6;     % max
kgrid = linspace(kmin,kmax,nbk)';

% initialization of various matrices
Tv = zeros(nbk,nbe);
dr = zeros(nbk,nbe);
C_fi = zeros(nbk,nbe);
util = zeros(nbk,nbe);

% computing the initial guess, by assuming deterministic shocks
temp = bsxfun(@times,muZ*(kgrid.^alpha) + (1-delta).*kgrid - kgrid,ones(nbk,nbe));
V = bsxfun(@rdivide, temp.^(1-sigma)-1, 1-sigma);
neg = temp<0;
V(neg) = -1e12;
%V = zeros(nbk,nbe);
% starting iterations
for it=1:3000
    
    for i=1:nbk    
        for k = 1:nbe
        C_fi (:,k) = Z(k)*(kgrid(i)^alpha) + (1-delta)*kgrid(i) - kgrid;      
        util (:,k) = (C_fi(:,k).^(1-sigma)-1)/(1-sigma);
        end
        neg = C_fi<0;
        util(neg) = -1e12;
        [Tv(i,:), dr(i,:)] = max (util + beta*(V*P));
    end
    
    crit = norm(Tv-V);
    if crit>epsi,V = Tv;
    else break;end
    
end

Kp_fi = kgrid(dr); % getting the optimal policy K'* for each K in the grid

% computing the optimal policy C* for each K in the grid
for k=1:nbe;
C_fi(:,k) = Z(k).*(kgrid.^alpha) + (1-delta).*kgrid - Kp_fi(:,k);
end
neg = C_fi<0;
C_fi(neg) = 1e-6;
util(neg) = -1e12;

t_fi = toc; % getting the duration in processor time

load Kp;
load C;
load t;

load Kp2;
load C2;
load t2;

norm (Kp_fi-Kp)
norm (C_fi-C)
t_fi-t

norm (Kp_fi-Kp2)
norm (C_fi-C2)
t_fi-t2
