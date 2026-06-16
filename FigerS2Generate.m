clear; clc; close all;
rng(20260516,'twister');   
%%
p.piS   = 10;
p.gamma = 0.04;
p.mu    = 0.01;
p.omega = 0.08;
p.tau   = 0.03;
p.muL   = 0.01;
p.piN   = 100;
p.alpha = 0.6;
p.ku    = 0.08;
p.eps   = 0.05;
p.ki    = 0.16;
p.muE0  = 0.01;
p.a     = 4.0;
p.beta0 = 0.45;
p.b     = 0.15;

xv = linspace(0,1,26);
rows = zeros(numel(xv)*20, 7);
rowi = 0;
%% inderect numerical solution to ODE system
for x = xv
    for rep = 1:20
        muE  = p.muE0*(1 + p.a*x)*lognmult(0.05);
        beta = p.beta0*(1 + p.b*x)*lognmult(0.05);

        y0 = [1000; 1; 0; 1000; 1];   % [S, I_E, I_L, N_u, N_i]
        opts = odeset('RelTol',1e-7,'AbsTol',1e-7);
        [t,y] = ode45(@(t,y) rhs_indirect(t,y,beta,muE,p), [0 600], y0, opts);

        S  = y(:,1);
        IE = y(:,2);
        IL = y(:,3);
        Nu = y(:,4);
        Ni = y(:,5);

        force = beta .* Ni ./ (Ni + Nu + 1e-12);

        incident = trapz(t, force .* S) * lognmult(0.03);
        reservoir_fraction = (Ni(end) / (Ni(end) + Nu(end) + 1e-12)) * lognmult(0.03);

        
        R0_like = (beta*p.alpha) / (muE + p.gamma + p.mu) * lognmult(0.03);

        rowi = rowi + 1;
        rows(rowi,:) = [x, rep, muE, beta, R0_like, incident, reservoir_fraction];
    end
end

df = array2table(rows, 'VariableNames', { ...
    'x','replicate','mu_E','beta','R0_like', ...
    'cumulative_incident_infections','final_contaminated_reservoir_fraction'});

writetable(df, 'indirect_reservoir_point_outputs.csv');
%%
% Figure
fig = figure('Color','w','Position',[100 100 1200 800]);
tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

pairs = { ...
    'beta', 'mu_E', 'Host infection rate, $\beta$', 'Observed early-stage virulence, $\mu_E$', 'A'; ...
    'R0_like', 'mu_E', 'Observed $R_0$', 'Observed early-stage virulence, $\mu_E$', 'B'; ...
    'beta', 'cumulative_incident_infections', 'Host infection rate, $\beta$', 'Observed cumulative infections', 'C'; ...
    'final_contaminated_reservoir_fraction', 'cumulative_incident_infections', 'Final contaminated-reservoir fraction', 'Observed cumulative infections', 'D' ...
};

cmap = [ ...
    0.00 0.45 0.74; ...
    0.85 0.33 0.10; ...
    0.47 0.67 0.19; ...
    0.49 0.18 0.56; ...
    0.30 0.30 0.30];

for k = 1:4
    ax = nexttile; hold(ax,'on');

    xcol = pairs{k,1};
    ycol = pairs{k,2};
    x = df.(xcol);
    y = df.(ycol);

    scatter(ax, x, y, 34, df.x, 'filled', ...
        'MarkerFaceAlpha',0.78, 'MarkerEdgeColor','none');
    % colormap(ax,cividis(256))
    % clim(ax,[0 1])
    idx = isfinite(x) & isfinite(y);
    z = polyfit(x(idx), y(idx), 1);
    xx = linspace(min(x), max(x), 100);
    plot(ax, xx, polyval(z,xx), 'Color',[0.25 0.25 0.25], 'LineWidth',1.8);

    r = corr(x(idx), y(idx));

    
text(ax,0.02,1.08,pairs{k,5}, ...
    'Units','normalized', ...
    'FontSize',22, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'Clipping','off');

    text(ax, 0.18, 1.04, sprintf('r = %.2f', r), ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontSize',20, ...
        'FontWeight','bold', ...
        'Color',[0.20 0.20 0.20], ...
        'Clipping','off');

    xlabel(ax, pairs{k,3}, 'Interpreter','latex', 'FontSize',18);
    ylabel(ax, pairs{k,4}, 'Interpreter','latex', 'FontSize',18);

    ax.FontSize = 18;
    ax.LineWidth = 1.0;
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.TickLength = [0.015 0.015];
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridAlpha = 0.18;
    ax.Colormap = cmap;
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Simulated variant index, $x$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 18;
cb.TickLabelInterpreter = 'latex';
%%
exportgraphics(fig, 'surface_model_figure.pdf', 'ContentType','vector');
exportgraphics(fig, 'surface_model_figure.png', 'Resolution',600);

fprintf('wrote indirect_reservoir_point_outputs.csv and surface_model_figure.pdf/png\n');
%% inderect ODE system
function dy = rhs_indirect(~, y, beta, muE, p)
S  = y(1);
IE = y(2);
IL = y(3);
Nu = y(4);
Ni = y(5);

%force = beta * Ni / (Ni + Nu + 1e-12);

dS  = p.piS + p.gamma*IE - beta*S*(Ni/(Ni+Nu+1e-12)) - p.mu*S;
dIE = beta*S*(Ni/(Ni+Nu+1e-12)) - (p.omega + p.tau + muE + p.gamma)*IE;
dIL = p.omega*IE - (p.muL + p.tau)*IL;
dNu = p.piN - p.alpha*(IE + IL)*Nu/(Ni + Nu + 1e-12) - p.ku*Nu + p.eps*Ni;
dNi = p.alpha*(IE + IL)*Nu/(Ni + Nu + 1e-12) - p.ki*Ni - p.eps*Ni;

dy = [dS; dIE; dIL; dNu; dNi];
end

function x = lognmult(cv)
sigma = sqrt(log(1 + cv^2));
mu = -0.5*sigma^2;
x = exp(mu + sigma*randn);
end