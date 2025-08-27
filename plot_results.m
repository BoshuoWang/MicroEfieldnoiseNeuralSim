load('UF_HH_3um\UF_HH_3um_compiled_result.mat');

UF_HH_3um = compiled_results;
UF_HH_3um.th_macro = reshape(UF_HH_3um.th_macro, [100, 7]);
UF_HH_3um.th_micro = reshape(UF_HH_3um.th_micro, [100, 7]);
UF_HH_3um.th_per_diff_micro = reshape(UF_HH_3um.th_per_diff_micro, [100, 7]);

load('UF_HH_0.3um\UF_HH_0.3um_compiled_result.mat');

UF_HH_03um = compiled_results;
UF_HH_03um.th_macro = reshape(UF_HH_03um.th_macro, [100, 7]);
UF_HH_03um.th_micro = reshape(UF_HH_03um.th_micro, [100, 7]);
UF_HH_03um.th_per_diff_micro = reshape(UF_HH_03um.th_per_diff_micro, [100, 7]);

load('PE_HH_3um\PE_HH_3um_compiled_result.mat');

PE_HH_3um = compiled_results;
PE_HH_3um.th_macro = reshape(PE_HH_3um.th_macro, [40, 7]);
PE_HH_3um.th_micro = reshape(PE_HH_3um.th_micro, [40, 7]);
PE_HH_3um.th_per_diff_micro = reshape(PE_HH_3um.th_per_diff_micro, [40, 7]);

load('PE_HH_0.3um\PE_HH_0.3um_compiled_result.mat');

PE_HH_03um = compiled_results;
PE_HH_03um.th_macro = reshape(PE_HH_03um.th_macro, [40, 7]);
PE_HH_03um.th_micro = reshape(PE_HH_03um.th_micro, [40, 7]);
PE_HH_03um.th_per_diff_micro = reshape(PE_HH_03um.th_per_diff_micro, [40, 7]);

for ii = size(UF_HH_3um.PW, 2) : -1 : 1
    if UF_HH_3um.PW(1, ii) < 1
        xticklabels{ii} = sprintf('$%g\\  \\mathrm{\\mu s}$', UF_HH_3um.PW(1, ii)*1000);
    else
        xticklabels{ii} = sprintf('$%g\\  \\mathrm{ms}$', UF_HH_3um.PW(1, ii));
    end
end
%%
load('UF_RMG_3um\UF_RMG_3um_compiled_result.mat');

UF_RMG_3um = compiled_results;
UF_RMG_3um.th_macro = reshape(UF_RMG_3um.th_macro, [100, 7]);
UF_RMG_3um.th_micro = reshape(UF_RMG_3um.th_micro, [100, 7]);
UF_RMG_3um.th_per_diff_micro = reshape(UF_RMG_3um.th_per_diff_micro, [100, 7]);

load('UF_RMG_0.3um\UF_RMG_0.3um_compiled_result.mat');

UF_RMG_03um = compiled_results;
UF_RMG_03um.th_macro = reshape(UF_RMG_03um.th_macro, [100, 7]);
UF_RMG_03um.th_micro = reshape(UF_RMG_03um.th_micro, [100, 7]);
UF_RMG_03um.th_per_diff_micro = reshape(UF_RMG_03um.th_per_diff_micro, [100, 7]);

load('PE_RMG_3um\PE_RMG_3um_compiled_result.mat');

PE_RMG_3um = compiled_results;
PE_RMG_3um.th_macro = reshape(PE_RMG_3um.th_macro, [40, 7]);
PE_RMG_3um.th_micro = reshape(PE_RMG_3um.th_micro, [40, 7]);
PE_RMG_3um.th_per_diff_micro = reshape(PE_RMG_3um.th_per_diff_micro, [40, 7]);

load('PE_RMG_0.3um\PE_RMG_0.3um_compiled_result.mat');

PE_RMG_03um = compiled_results;
PE_RMG_03um.th_macro = reshape(PE_RMG_03um.th_macro, [40, 7]);
PE_RMG_03um.th_micro = reshape(PE_RMG_03um.th_micro, [40, 7]);
PE_RMG_03um.th_per_diff_micro = reshape(PE_RMG_03um.th_per_diff_micro, [40, 7]);

%%
col = lines(7);

clearvars h*
for hh = 1 : 2
    h_f(hh) = figure('Color', 'w', 'Position', [75 + 400 * hh, 50, 640, 480]);
    h_a(hh) = axes(h_f(hh), 'Tag', 'Main plot');
    h_ax_1(hh) = axes(h_f(hh), 'Position',[0.005, 0.945, 0.05, 0.05], 'Clipping', 'off', 'Box', 'off', 'Tag', 'Figure label');
    axis(h_ax_1(hh), 'off');
    axis(h_ax_1(hh), [0, 1, -1, 0]);
end
text(h_ax_1(1), 0, 0, 0, '$\mathbf{B_{iv}}$', 'Interpreter', 'Latex', 'FontSize', 30, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
text(h_ax_1(2), 0, 0, 0, '$\mathbf{C_{iv}}$', 'Interpreter', 'Latex', 'FontSize', 30, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');

set(h_a, 'NextPlot', 'Add', 'TickLabelInterpreter', 'Latex', 'FontSize', 16, ...
    'Box', 'On', 'LineWidth', 1, 'TickDir', 'both', 'TickLength', [0.005, 0.01]);

% h_b1{1} = boxplot(h_a(1), UF_HH_3um.th_per_diff_micro, 'plotstyle', 'compact', 'Positions', (1:7)-0.1, 'Colors', col(1, :), 'jitter', 0.25);
% h_b1{2} = boxplot(h_a(1), UF_HH_03um.th_per_diff_micro, 'plotstyle', 'compact', 'Positions', (1:7)+0.1, 'Colors', col(2, :), 'jitter', 0.25);
% h_b1{3} = boxplot(h_a(1), UF_HH_30nm.th_per_diff_micro, 'plotstyle', 'compact', 'Positions', (1:7)+0.25, 'Colors', col(3, :), 'jitter', 0.1);

h_b1{1} = boxplot(h_a(1), UF_HH_3um.th_per_diff_micro, 'plotstyle', 'compact', 'Positions', (1:7)-0.15, 'Colors', col(1, :), 'jitter', 0.25);
h_b1{2} = boxplot(h_a(1), UF_HH_03um.th_per_diff_micro, 'plotstyle', 'compact', 'Positions', (1:7)-0.05, 'Colors', col(6, :), 'jitter', 0.25);
h_b1{3} = boxplot(h_a(1), UF_RMG_3um.th_per_diff_micro, 'plotstyle', 'compact', 'Positions', (1:7)+0.05, 'Colors', col(7, :), 'jitter', 0.25);
h_b1{4} = boxplot(h_a(1), UF_RMG_03um.th_per_diff_micro, 'plotstyle', 'compact', 'Positions', (1:7)+0.15, 'Colors', col(2, :), 'jitter', 0.25);
% set(handle(h_b), 'LineWidth', 1)

xlabel(h_a(1), {'Pulse duration'}, 'FontSize', 18, 'Interpreter', 'latex');
ylabel(h_a(1), 'Relative threshold difference', 'FontSize', 18, 'Interpreter', 'latex');
% ytick = h_a(1).YTick;
ylim(h_a(1), [-3.75, 3.75]);
ytick = -5:1:5;
yticklabels = cell(size(ytick));
for ii = 1 : length(ytick)
    yticklabels{ii} = sprintf('$$%g\\%%$$', ytick(ii));
end
set(h_a(1), 'XTick', 1:7, 'XTickLabel', xticklabels, 'YTickLabel', yticklabels, 'YTick', ytick, 'FontSize', 14, 'TickLabelInterpreter', 'Latex');
h_line(1) = plot(h_a(1), xlim(h_a(1)), 0.05*[1,1], '--', 'Color', [1,1,1]*0.5, 'LineWidth', 1);
h_line(2) = plot(h_a(1), xlim(h_a(1)), -0.05*[1,1], '--', 'Color', [1,1,1]*0.5, 'LineWidth', 1);
uistack(h_line, "bottom");
ylim(h_a(1), [-4, 4]);
h_line2(1) = plot(h_a(2), [-1, -1], [0, 1], '-', 'LineWidth', 2.5, 'color', col(1, :));
h_line2(2) = plot(h_a(2), [-1, -1], [0, 1], '-', 'LineWidth', 2.5, 'color', col(6, :));
h_line2(3) = plot(h_a(2), [-1, -1], [0, 1], '-', 'LineWidth', 2.5, 'color', col(7, :));
h_line2(4) = plot(h_a(2), [-1, -1], [0, 1], '-', 'LineWidth', 2.5, 'color', col(2, :));


h_b2{1} = boxplot(h_a(2), PE_HH_3um.th_per_diff_micro, 'plotstyle', 'compact', 'Positions', (1:7)-0.15, 'Colors', col(1, :), 'jitter', 0.25);
h_b2{2} = boxplot(h_a(2), PE_HH_03um.th_per_diff_micro, 'plotstyle', 'compact', 'Positions', (1:7)-0.05, 'Colors', col(6, :), 'jitter', 0.25);
h_b2{3} = boxplot(h_a(2), PE_RMG_3um.th_per_diff_micro, 'plotstyle', 'compact', 'Positions', (1:7)+0.05, 'Colors', col(7, :), 'jitter', 0.25);
h_b2{4} = boxplot(h_a(2), PE_RMG_03um.th_per_diff_micro, 'plotstyle', 'compact', 'Positions', (1:7)+0.15, 'Colors', col(2, :), 'jitter', 0.25);



xlabel(h_a(2), {'Pulse duration'}, 'FontSize', 18, 'Interpreter', 'latex');
ylabel(h_a(2), 'Relative threshold difference', 'FontSize', 18, 'Interpreter', 'latex');
ylim(h_a(2), [-4, 4]);
ytick = -4:1:4;
yticklabels = cell(size(ytick));
for ii = 1 : length(ytick)
    yticklabels{ii} = sprintf('$$%g\\%%$$', ytick(ii));
end
set(h_a(2), 'XTick', 1:7, 'XTickLabel', xticklabels, 'YTickLabel', yticklabels, 'YTick', ytick, 'FontSize', 14, 'TickLabelInterpreter', 'Latex');
h_line(1) = plot(h_a(2), xlim(h_a(2)), 0.05*[1,1], '--', 'Color', [1,1,1]*0.5, 'LineWidth', 1);
h_line(2) = plot(h_a(2), xlim(h_a(2)), -0.05*[1,1], '--', 'Color', [1,1,1]*0.5, 'LineWidth', 1);
uistack(h_line, "bottom");


xlim(h_a, [0.5, 7.5]);
h_a(1).Position = [0.17, 0.15, 0.78, 0.78];
h_a(2).Position = [0.17, 0.15, 0.78, 0.78];


h_b =  handle([h_b1{1},h_b1{2}, h_b2{1},h_b2{2}]);
set(findobj(h_b, 'Tag', 'Whisker'),       'LineWidth', 1.5);
set(findobj(h_b, 'Tag', 'Box'),           'LineWidth', 5);
set(findobj(h_b, 'Tag', 'MedianOuter'),   'Marker', 'o', 'MarkerSize', 7, 'LineWidth', 1);
set(findobj(h_b, 'Tag', 'Outliers'),      'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'None');

h_a(1).XLabel.Units = 'normalized';
h_a(1).XLabel.Position(2) = h_a(1).XLabel.Position(2) - 0.05;

h_a(2).XLabel.Units = 'normalized';
h_a(2).XLabel.Position(2) = h_a(2).XLabel.Position(2) - 0.05;

%%
h_ax_in = axes(h_f(1), 'Position', [0.5, 0.892, 0.35, 0.0960], 'Color', 'w', 'NextPlot', 'Add');

cla(h_ax_in);
set(h_ax_in, 'XColor', 'none', 'YColor', 'none',  'Position', [0.5, 0.75, 0.4, 0.15], 'Tag', 'inset');
axis(h_ax_in,[-0.05,0.95,0.05,0.85]);

text(h_ax_in, -0.1, 0.5, 'Unmyelinated',  'Interpreter', 'Latex', 'HorizontalALignment', 'Left');
text(h_ax_in, -0.1, 0.2, 'Myelinated',  'Interpreter', 'Latex', 'HorizontalALignment', 'Left');
text(h_ax_in, 0.5, 0.7, '$3 \ \mathrm{\mu m}$',  'Interpreter', 'Latex', 'HorizontalALignment', 'Center');
text(h_ax_in, 0.8, 0.7, '$0.3 \ \mathrm{\mu m}$',  'Interpreter', 'Latex', 'HorizontalALignment', 'Center');
set(findobj(h_ax_in,'Type', 'text'), 'VerticalAlignment', 'Middle', 'FontSize', 14, 'Color', 'k');


plot(h_ax_in, [0.3,0.5]+0.1, [0.5,0.5], '-', 'LineWidth', 1.5, 'Color', col(1, :));
plot(h_ax_in, [0.6,0.8]+0.1, [0.5,0.5], '-', 'LineWidth', 1.5, 'Color', col(6, :));
plot(h_ax_in, [0.3,0.5]+0.1, [0.2,0.2], '-', 'LineWidth', 1.5, 'Color', col(7, :));
plot(h_ax_in, [0.6,0.8]+0.1, [0.2,0.2], '-', 'LineWidth', 1.5, 'Color', col(2, :));
uistack(h_f(1).Children(1), 'down')
% h_f(1).Children
%%

cla(h_ax_1(1));
cla(h_ax_1(2));
text(h_ax_1(1), 0, 0, 0, '$\mathbf{A}$', 'Interpreter', 'Latex', 'FontSize', 30, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
text(h_ax_1(2), 0, 0, 0, '$\mathbf{B}$', 'Interpreter', 'Latex', 'FontSize', 30, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
saveas(h_f(1), fullfile('Figures_Full', '3A.fig'));
saveas(h_f(2), fullfile('Figures_Full', '3B.fig'));

%%
% set(findobj(h_b, 'Tag', 'Whisker'),'LineWidth', 2);
%         set(findobj(h_b, 'Tag', 'Box'),'LineWidth', 7);
%         set(findobj(h_b, 'Tag', 'Outliers'), 'Marker', 'none');