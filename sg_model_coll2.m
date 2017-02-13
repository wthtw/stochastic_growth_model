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
nbe = 15;           % number of shocks
se = 0.12;          % volatility of the shock
mu = 0;             % mean of shocks residuals iid normal process

% Matrice of transition probabilities using tauchen and hussey
muZ = 0.5;            % mean of shocks
[Z,P] = tauchenhussey(nbe,muZ,rho,se,se);

% constructing grid values for K, using collocation method
nbk     = 50;        % number of collocation points
kmin  = 0.2;
kmax  = 6;
basis = fundefn('lin',nbk,kmin,kmax);  % cubic spline basis   
kgrid = funnode(basis);                    
Phi   = funbas(basis);

% setting initial guess for basis coefficients c
c = zeros(nbk,nbe);                     
%load c_gscoll;

% initialization of various matrices
Tv = zeros(nbk,nbe); 
dr = zeros(nbk,nbe);
C = zeros(nbk,nbe);
util = zeros(nbk,nbe);

for it=1:4000
    
cold   = c;

for p=1:nbk
temp1=(kgrid(p).^alpha)*(Z') + (1-delta).*kgrid(p);
C = bsxfun(@plus,temp1,-kgrid);
util = bsxfun(@rdivide, C.^(1-sigma)-1, 1-sigma); % bsxfun do element-by-element operations between arrays, here a division
neg = C<0;
util(neg) = -1e12;
temp2 = bsxfun(@plus, util, beta*(funeval(c,basis,kgrid)*P));
[Tv(p,:), dr(p,:)] = max(temp2);
end

c = Phi\Tv;
crit = norm(c-cold);
if crit<epsi, break,end

end

Kp2 = kgrid(dr);
temp2=bsxfun(@times,(kgrid.^alpha),Z');
temp3=bsxfun(@plus,temp2,(1-delta).*kgrid);
C2 = bsxfun(@plus,temp3,-Kp2);
neg = C2<0;
C2(neg) = 1e-6;
util(neg) = -1e12;

t2=toc;

save ('Kp2','Kp2');
save ('C2','C2');
save ('t2','t2');
save ('c','c');
