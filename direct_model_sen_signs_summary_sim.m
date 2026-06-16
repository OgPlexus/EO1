
clear; clc;close all

df = readtable('direct_model_sensitivity.csv');

relationships = { ...
    'r_mus_beta', '$\mu_s$ vs $\beta$'; ...
    'r_mus_R0',   '$\mu_s$ vs $R_0$'; ...
    'r_D_beta',   '$D(\infty)$ vs $\beta$'; ...
    'r_D_final',  '$D(\infty)$ vs final epidemic size'};
%%

rel = strings(4,1);
negative = zeros(4,1);
near_zero = zeros(4,1);
positive = zeros(4,1);
n = zeros(4,1);

for i = 1:4
    col = relationships{i,1};
    label = relationships{i,2};
    vals = df.(col);
    vals = vals(~isnan(vals));

    rel(i) = label;
    negative(i) = sum(vals < -0.2);
    near_zero(i) = sum(vals >= -0.2 & vals <= 0.2);
    positive(i) = sum(vals > 0.2);
    n(i) = numel(vals);
end
%%
out = table(rel, negative, near_zero, positive, n, ...
    'VariableNames', {'relationship','negative','near_zero','positive','n'});

writetable(out, 'direct_model_sensitivity_signs_summary.csv');
