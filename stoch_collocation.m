%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   stoch_collocation.m:  A Matlab program to solve a simple 
%   stochastic growth model via collocation.
%
%   Youssef de Madeen Amadou, Winter 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializing time chrono
tic;
clc;clear;

%
%  Parameter values
%

sigma = 1.50;       % utility parameter
delta = 0.10;       % depreciation rate
beta = 0.95;        % discount factor
alpha = 0.30;       % capital elasticity of output
rho = 0.80;         % persistence of the shock on Z
epsi = 1e-5;        % convergence parameter

% Discretization of the shocks, Markov AR(1) process
nbe = 15;           % number of shocks, 15
se = 0.12;          % volatility of the shock
mu = 0;             % mean of shocks residuals iid normal process

% Matrice of transition probabilities using tauchen and hussey
muZ = 0.5;            % mean of shocks
[Z,P] = tauchenhussey(nbe,muZ,rho,se,se);

% Constructing grid values for K, using collocation method
nbk     = 50;        % number of collocation points, 50
kmin  = 0.2;
kmax  = 6;
basis = fundefn('lin',nbk,kmin,kmax);  % cubic spline basis   
kgrid = funnode(basis);                    
Phi   = funbas(basis);

% Setting initial guess for basis coefficients c
c = zeros(nbk,nbe);
%load c;
%load c_gscoll;

% Initializating some matrices
Tv = zeros(nbk,nbe); 
dr = zeros(nbk,nbe);
C = zeros(nbk,nbe);
util = zeros(nbk,nbe);

% Iterations
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

% Getting the optimal policy K'* for each K in the grid
Kp = kgrid(dr); 

% Computing the optimal policy C* for each K in the grid
for k=1:nbe;
C(:,k) = Z(k).*(kgrid.^alpha) + (1-delta).*kgrid - Kp(:,k);
end
neg = C<0;
C(neg) = 1e-6;
util(neg) = -1e12;

% Closing chrono
t = toc;

% Saving results
save ('Kp','Kp'); save ('C','C'); save ('t','t'); save ('c','c')


%% Simulating the solved stochastic growth model

T = 200;        % number of periods
p = rand(T,1);  % uniformly distributed 'p' considered as probabilities
Zsim = zeros(T+1,1);
Zsim(1) = Z(9); % starting value for Z

% Simulating the shocks, using probabilities in p
for j=2:T+1
        l = find(Z==Zsim(j-1));
        m = find(p(j-1)<=cumsum(P(l,:),2),1,'first');
        if find(p(j-1)>cumsum(P(l,:),2),1,'last') == m - 1;
            Zsim(j)=Z(m);
        end
end

Ksim = zeros(T+1,1);
Ksim(1) = kgrid(kmax);
Ysim = zeros(T,1);

% Simulating Stock of capital and Production evolution
% Knowing Zt-1, we interpolate the appropriate decision rule fonction for
% Kt on the grid
for t=1:T;
Ysim(t) = Zsim(t)*(Ksim(t).^alpha);
l = find(Z==Zsim(t));
Ksim(t+1) = interp1(kgrid,Kp2(:,l),Ksim(t),'spline');  
end

% Computing implied Investment and Consumption
Isim = Ksim(2:T+1)-(1-delta)*Ksim(1:T);
Csim = Ysim-Isim;

figure;
subplot(2,2,1), plot(Csim,'r'), title('Consumption'), xlabel('Time');
subplot(2,2,2), plot(Ysim,'b'), title('Production'), xlabel('Time');
subplot(2,2,3), plot(Ksim,'r'), title('Stock of Capital'), xlabel('Time'), axis([1 T+1 0.2 6]);
subplot(2,2,4), plot(Isim,'b'), title('Investment'), xlabel('Time');
