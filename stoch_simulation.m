%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   stoch_simulation.m:  A Matlab program to simulate a simple stochastic 
%   growth model using Function Iteration
%
%   Youssef de Madeen Amadou, Winter 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
