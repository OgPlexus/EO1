

clear; close all; clc; 

pts = readtable(fullfile('figure1_theoretical_projection.csv')); % toy data for E-O is in this file
%%
panels = {'A','B','C','D'};
maps   = {'$V=E,\ T=O$','$V=O,\ T=O$','$V=E,\ T=E$','$V=O,\ T=E$'};
vsrc   = {'E','O','E','O'};
tsrc   = {'O','O','E','E'};

out = table();
for k = 1:numel(panels)
    n = height(pts);
    tmp = table(repmat(string(panels{k}),n,1), repmat(string(maps{k}),n,1), ...
        pts.genotype_id, pts.E, pts.O, pts.(vsrc{k}), pts.(tsrc{k}), ...
        'VariableNames', {'panel','mapping','genotype_id','E','O','V','T'});
    out = [out; tmp]; 
end
%%

writetable(out, fullfile('figure1_projected_observables.csv'));
%%

fig = figure('Color','w','Position',[100 100 1080 770]);
%sgtitle('A single efficiency-opportunity constraint can produce different observed virulence-transmission patterns.','FontSize',13);

ax0 = subplot(2,3,[1 4]);
ax0.Position = [0.045 0.13 0.3 0.65];
set(ax0,'Color',[0.96 0.96 0.96]); hold(ax0,'on');
plot(ax0, pts.E, pts.O, 'k-', 'LineWidth', 2.8);
scatter(ax0, pts.E, pts.O, 58, 'k', 'filled');
axis(ax0,[0 1 0 1]); axis(ax0,'square');
xlabel(ax0,'Efficiency, $E$','FontSize',11.5, 'Interpreter','latex');
ylabel(ax0,'Opportunity, $O$','FontSize',11.5,'Interpreter','latex');
title(ax0,'Latent $E-O$','FontSize',12.5,'Interpreter','latex');
set(ax0,'FontSize',14,'LineWidth',1.25);
text(ax0,-0.05,1.05,'A','Units','normalized','FontSize',14,'FontWeight','bold');

tilepos = [2 3 5 6];
for i = 1:4
    ax = subplot(2,3,tilepos(i));
    sub = out(out.panel == string(panels{i}),:);
    r = corr(sub.V, sub.T);
    [~,idx] = sort(sub.V);
    hold(ax,'on');
    plot(ax, sub.V(idx), sub.T(idx), '-', 'Color',[0.02 0.02 0.02], 'LineWidth',1.05);
    scatter(ax, sub.V, sub.T, 30, [0.01 0.01 0.01], 'filled', 'MarkerFaceAlpha',0.62);
    axis(ax,[0 1 0 1]); axis(ax,'square');
    xlabel(ax,'Virulence, $V$','FontSize',10,'Interpreter','latex');
    ylabel(ax,'Transmission, $T$','FontSize',10,'Interpreter','latex');
    title(ax,maps{i},'Interpreter','latex','FontSize',11.5,'Color',[0.15 0.15 0.15]);
    set(ax,'FontSize',14,'LineWidth',0.75,'XColor',[0.30 0.30 0.30],'YColor',[0.30 0.30 0.30]);
    text(ax,-0.05,1.07,char('B'+(i-1)),'Units','normalized','FontSize',14,'FontWeight','bold', ...
        'Color',[0.10 0.10 0.10]);
idxMark = round(0.85*numel(idx));  
x0 = sub.V(idx(idxMark));
y0 = sub.T(idx(idxMark));

dx = 0.09;
dy = 0.02;

scatter(ax, x0+dx, y0+dy, 260, 'o', 'MarkerEdgeColor',[0.15 0.15 0.15],'MarkerFaceColor','none','LineWidth',1.1);

text(ax, x0+dx, y0+dy, ternary(r>0,'$+$','$-$'), ...
    'Interpreter','latex', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',18, ...
    'FontWeight','bold', ...
    'Color',[0.10 0.10 0.10]);
end
%%
set(gcf,'PaperPositionMode','auto');
exportgraphics(fig, fullfile('figure1_projections.pdf'), 'ContentType','vector');
exportgraphics(fig, fullfile('figure1_projections.png'), 'Resolution',300);
close(fig);


function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
end
