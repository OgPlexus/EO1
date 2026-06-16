clear all; close all; clc;
rng(20260516,'twister');

out = fullfile('direct_transmission_point_outputs.csv');
%%

mu = 0; eps = 1/3; omega = 1/4; p = 0.35;
nu0 = 1/7; mus0 = 0.004; alpha0 = 0.12; beta0 = 0.45;
a = 5.0; bA = 0.8; bB = 1.4; d = 4.0; N = 100000;

xv = linspace(0,1,26);
n = numel(xv)*20;

X = zeros(n,1); Rep = zeros(n,1); Alpha = zeros(n,1); Beta = zeros(n,1);
Nu = zeros(n,1); MuS = zeros(n,1); R0 = zeros(n,1); Deaths = zeros(n,1);
FinalSize = zeros(n,1); PeakInf = zeros(n,1);

k = 0;
opts = odeset('RelTol',1e-7,'AbsTol',1e-7);

for x = xv
    for rep = 1:20
        k = k + 1;

        alpha = alpha0*(1 + bA*x)*lognormcv(0.05);
        beta  = beta0 *(1 + bB*x)*lognormcv(0.05);
        nu    = nu0   *(1 + d*x)*lognormcv(0.05);
        mus   = mus0  *(1 + a*x)*lognormcv(0.05);

        y0 = [N-1; 0; 0; 1; 0; 0];
        [~, y] = ode45(@(t,y) rhs(t,y,alpha,beta,eps,omega,p,nu,mus,mu,N), [0 500], y0, opts);

        X(k) = x; Rep(k) = rep; Alpha(k) = alpha; Beta(k) = beta;
        Nu(k) = nu; MuS(k) = mus;
        FinalSize(k) = (y(end,5) + y(end,6))*lognormcv(0.03);
        Deaths(k) = y(end,6)*lognormcv(0.03);
        PeakInf(k) = max(y(:,3) + y(:,4))*lognormcv(0.03);
        R0(k) = ((alpha + beta)/(nu + mus + mu))*lognormcv(0.03);
    end
end

T = table(X,Rep,Alpha,Beta,Nu,MuS,R0,Deaths,FinalSize,PeakInf, ...
    'VariableNames', {'x','replicate','alpha','beta','nu','mu_s', ...
    'R0_summary','cumulative_deaths','final_size','peak_infectious'});

writetable(T,out); % save the data
%%
function r = lognormcv(cv)
s = sqrt(log(1 + cv^2));
r = exp(-0.5*s^2 + s*randn);
end
%direct SEIIR ODE system
function dydt = rhs(~,y,alpha,beta,eps,omega,p,nu,mus,mu,N)
S = y(1); E = y(2); IA = y(3); IS = y(4); R = y(5); D = y(6);
lam = (alpha*IA + beta*IS)/N;

dydt = [ mu*(N-S) - lam*S
         lam*S - (eps+mu)*E
         eps*E - (omega+mu)*IA
         (1-p)*omega*IA - (nu+mus+mu)*IS
         p*omega*IA + nu*IS - mu*R
         mus*IS ];
end