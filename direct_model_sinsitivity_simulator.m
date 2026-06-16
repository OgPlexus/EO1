
clear; clc;close all;
rng(20260516,'twister');
%%
mu    = 0.01;
eps   = 1/3;
omega = 1/4;
p     = 0.35;
nu0   = 1/7;
mus0  = 0.004;
alpha0= 0.12;
beta0 = 0.45;
N     = 100000;

xv   = linspace(0,1,16);
nRun = 150;
out  = zeros(nRun, 8);

opts = odeset('RelTol',1e-7,'AbsTol',1e-7);

for run = 1:nRun
    a  = 1 + 7*rand;       % U(1,8)
    bA = 2*rand;           % U(0,2)
    bB = 0.2 + 2.3*rand;   % U(0.2,2.5)
    d  = 7*rand;           % U(0,7)

    mus   = mus0   * (1 + a*xv);
    beta  = beta0  * (1 + bB*xv);
    alpha = alpha0 * (1 + bA*xv);
    nu    = nu0    * (1 + d*xv);

    R0 = (alpha + beta) ./ (nu + mus + mu);   

    D = zeros(size(xv));
    F = zeros(size(xv));   

    for j = 1:numel(xv)
        y0 = [N-1; 0; 0; 1; 0; 0];  
        [~, y] = ode45(@(t,y) rhs_direct(t,y,alpha(j),beta(j),eps,omega,p,nu(j),mus(j),mu,N), ...
                       [0 500], y0, opts);

        D(j) = y(end,6);
        F(j) = y(end,5) + y(end,6);
    end

    out(run,:) = [ ...
        a, bA, bB, d, ...
        pearsonr(mus, beta), ...
        pearsonr(mus, R0), ...
        pearsonr(D, beta), ...
        pearsonr(D, F)];
end
%%
T = array2table(out, 'VariableNames', ...
    {'a','bA','bB','d','r_mus_beta','r_mus_R0','r_D_beta','r_D_final'});

writetable(T, 'direct_model_sensitivity.csv');

%% ODE System
function dy = rhs_direct(~, y, alpha, beta, eps, omega, p, nu, mus, mu, N)
S  = y(1);
E  = y(2);
IA = y(3);
IS = y(4);
R  = y(5);
D  = y(6);

lam = (alpha*IA + beta*IS) / N;

dy = zeros(6,1);
dy(1) = mu*(N - S) - lam*S;
dy(2) = lam*S - (eps + mu)*E;
dy(3) = eps*E - (omega + mu)*IA;
dy(4) = (1-p)*omega*IA - (nu + mus + mu)*IS;
dy(5) = p*omega*IA + nu*IS - mu*R;
dy(6) = mus*IS;
end
%% corr r value %%%
function r = pearsonr(x, y)
c = corrcoef(x(:), y(:));
r = c(1,2);
end