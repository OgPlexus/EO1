
clear; clc; close all

sens = readtable('direct_model_sensitivity_signs_summary.csv');
%%
neg  = sens.negative;
near = sens.near_zero;
pos  = sens.positive;
y    = 1:height(sens);

fig = figure('Color','w','Position',[100 100 1200 680]);
ax = axes(fig); hold(ax,'on');

% Colorblind-friendly palette
cNeg  = [0.12 0.47 0.71];   % blue
cNear = [0.78 0.78 0.78];   % light gray
cPos  = [0.95 0.62 0.12];   % orange

b = barh(ax, y, [neg near pos], 'stacked', 'BarWidth', 0.55);
b(1).FaceColor = cNeg;  b(1).EdgeColor = 'k'; b(1).LineWidth = 0.7;
b(2).FaceColor = cNear; b(2).EdgeColor = 'k'; b(2).LineWidth = 0.7;
b(3).FaceColor = cPos;  b(3).EdgeColor = 'k'; b(3).LineWidth = 0.7;

for i = 1:height(sens)
    vals   = [neg(i), near(i), pos(i)];
    starts = [0, neg(i), neg(i) + near(i)];

    for j = 1:3
        if vals(j) > 0
            x = starts(j) + vals(j)/2;

            txtColor = 'k';
            if j == 1 && vals(j) >= 25
                txtColor = 'w';
            end

            text(ax, x, y(i), sprintf('%d', vals(j)), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontSize',14, ...
                'FontWeight','bold', ...
                'Color',txtColor);
        end
    end
end

ax.YTick = y;
ax.YTickLabel = sens.relationship;
ax.XLim = [0 max(sens.n)];
ax.YDir = 'reverse';
ax.Box = 'off';
ax.LineWidth = 1.0;
ax.FontSize = 15;
ax.FontName = 'Helvetica';
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0 0];
ax.XGrid = 'on';
ax.GridAlpha = 0.18;
ax.GridLineStyle = '-';

xlabel(ax, 'Number of sampled parameter sets', ...
    'Interpreter','latex', 'FontSize',18);

lgd = legend(ax, b, { ...
    'Negative ($r<-0.2$)', ...
    'Near zero ($-0.2 \leq r \leq 0.2$)', ...
    'Positive ($r>0.2$)'}, ...
    'Interpreter','latex', ...
    'Orientation','horizontal', ...
    'Box','off');

lgd.FontSize = 15;
lgd.Units = 'normalized';
lgd.Position = [0.24 0.86 0.52 0.08];
%%
set(fig, 'PaperPositionMode', 'auto');
exportgraphics(fig, 'direct_model_sensitivity_signs.pdf', 'ContentType','vector');
exportgraphics(fig, 'direct_model_sensitivity_signs.png', 'Resolution',300);