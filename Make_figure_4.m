clear all;close all;clc;
df = readtable('direct_transmission_point_outputs.csv');

pairs = {
    'B','beta','mu_s','Symptomatic transmission, $\beta$','Disease mortality, $\mu_s$','Mortality vs. transmission rate'
    'C','R0_summary','mu_s','$R_0$-like summary','Disease mortality, $\mu_s$','Mortality vs. R_0'
    'D','beta','cumulative_deaths','Symptomatic transmission, $\beta$','Cumulative deaths','Deaths vs. transmission rate'
    'E','final_size','cumulative_deaths','Final epidemic size','Cumulative deaths','Deaths vs. final epidemic size'
};
colors = [0.20 0.45 0.70;   % B
          0.85 0.33 0.10;   % C
          0.47 0.67 0.19;   % D
          0.49 0.18 0.56];  % E
%%
close all;
r = zeros(4,1);
for i = 1:4, r(i) = corr(df.(pairs{i,2}), df.(pairs{i,3})); end

f = figure('Color','w','Position',[100 100 1350 920]);
t = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

ax = nexttile([1 2]);
for i = 1:4
    barh(ax, i, r(i), 'FaceColor', colors(i,:), 'EdgeColor', 'none'); hold(ax,'on')
end
xline(ax,0,'k','LineWidth',1.1)
for i = 1:4
    text(ax,-1.05,i,pairs{i,1},'FontWeight','bold','FontSize',18,'HorizontalAlignment','left')
    if r(i) < 0
        hal = 'right'; dx = -0.04;
    else
        hal = 'left'; dx = 0.04;
    end
    text(ax,r(i)+dx,i,sprintf('%.2f',r(i)), ...
        'FontWeight','bold','FontSize',18,'HorizontalAlignment',hal)
end
set(ax,'YDir','reverse','YTick',[],'XLim',[-1.1 1.1],'FontName','Arial','FontSize',11,'LineWidth',1,'Box','off')
grid(ax,'on'); ax.GridAlpha = 0.15; ax.XColor = [0.2 0.2 0.2];
xlabel(ax,'Correlation across simulated epidemics, r','FontSize',18)
title(ax,'A. Correlations for four endpoint pairings','FontWeight','bold','FontSize',18,'HorizontalAlignment','center')%,'Interpreter','latex')

for i = 1:4
    ax = nexttile;
    x = df.(pairs{i,2}); y = df.(pairs{i,3});
    scatter(ax, x, y, 18, 'filled', ...
        'MarkerFaceColor', colors(i,:), 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none'); hold(ax,'on')
    p = polyfit(x, y, 1); xx = linspace(min(x), max(x), 100);
    plot(ax, xx, polyval(p, xx), 'Color', colors(i,:), 'LineWidth', 2.2)
    set(ax,'FontName','Arial','FontSize',18,'LineWidth',1,'Box','off')
    grid(ax,'on'); ax.GridAlpha = 0.15;
    title(ax,sprintf('%s. %s (r = %.2f)',pairs{i,1},pairs{i,6},r(i)), ...
        'FontWeight','bold','FontSize',18,'HorizontalAlignment','center')%,'Interpreter','latex')
    xlabel(ax,pairs{i,4},'FontSize',18,'Interpreter','latex')
    ylabel(ax,pairs{i,5},'FontSize',18,'Interpreter','latex')
end
%%
%exportgraphics(f,'adapted_model_figure.png','Resolution',1200)
exportgraphics(f,'adapted_model_figure.pdf','ContentType','vector')
