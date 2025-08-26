load("pSE_distribution.mat");
opt = optimset(optimset('fminsearch'), 'Display', 'on', ...
    'MaxFunEvals', 100000, 'MaxIter', 20000, ...
    'TolX', 1e-10, 'TolFun', 1e-8); % Options for fminsearch

dz = 0.1e-4;        % 1 um in cm
z_end = 1;          % 1 cm  
z_vec = 0:dz:z_end;

N = 10;

N_lambda = 25;
lambda_min = 1e-4;
lambda_mid = 15e-4;
lambda_max = 50e-4;
E_noise_max = 0.3;

dSE = mean(diff(SE_vec));
SE_vec2 = [(-1 : dSE : -dSE)' ; SE_vec];
SEinterp_vec2 = (SE_vec2(1:end-1) + SE_vec2(2:end))/2;


lambda = cell(N, 1);
mean_x = cell(N, 1);
pSE_actual = cell(N, 1);

for nn = N : -1 : 1
    T_start = tic;
    lambda{nn} = [rand(1, N_lambda-5) * (lambda_mid - lambda_min) + lambda_min, rand(1, 5) * (lambda_max - lambda_mid) + lambda_mid];
    phi{nn} = 2*pi * rand(size(lambda{nn}));
    
    
    A{nn} = rand(size(lambda{nn})) * E_noise_max;
    params_init = sort(A{nn});
    
    [params, RSS{nn}, exit_flag] = ...
        fminsearch(@(x) RSSerror(x, N_lambda, lambda{nn}, phi{nn}, z_vec, SE_vec, pSEinterp_vec), params_init, opt);
    
    RSS{nn} = sqrt(RSS{nn}*mean(diff(SEinterp_vec)));
    
    A{nn} = sort(abs(params));
    
    E_z{nn} = ones(size(z_vec));        % mV/cm
       
    for ii = 1 : N_lambda
        E_z{nn} = E_z{nn} + A{nn}(ii) * sin( 2*pi * z_vec/lambda{nn}(ii) - phi{nn}(ii));
    end
    
    pSE_actual{nn} = histcounts(abs(E_z{nn}), SE_vec, 'Normalization', 'pdf')';
    pSE_actual_signed{nn} = histcounts(E_z{nn}, SE_vec2, 'Normalization', 'pdf')';
    
    mean_x{nn} = trapz(SEinterp_vec, pSE_actual{nn} .* SEinterp_vec);
    T = toc(T_start);
    fprintf("Loop %d. Time: %g s. Exitflag: %d; fmin RSS: %g; Mean: %g\n", nn, T, exit_flag, RSS{nn}, mean_x{nn});
end
%%
pSE_actual_mean = mean([pSE_actual{:}], 2);
pSE_actual_std = std([pSE_actual{:}], 0, 2);

pSE_actual_signed_mean = mean([pSE_actual_signed{:}], 2);
pSE_actual_signed_std = std([pSE_actual_signed{:}], 0, 2);

    
clearvars h*
for hh = 1 : 2
    h_f(hh) = figure('Color', 'w', 'Position', [75 * hh, 50, 825, 720]);
    h_a(hh) = axes(h_f(hh));
    h_ax_1(hh) = axes(h_f(hh), 'Position',[0.005, 0.945, 0.05, 0.05],'Clipping','off','Box','off');
    axis(h_ax_1(hh), 'off');
    axis(h_ax_1(hh), [0, 1, -1, 0]);
end
text(h_ax_1(1), 0, 0, 0, '$\mathbf{A_{ii}}$', 'Interpreter', 'Latex', 'FontSize', 36, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
text(h_ax_1(2), 0, 0, 0, '$\mathbf{A_{iii}}$', 'Interpreter', 'Latex', 'FontSize', 36, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');

set(h_a, 'NextPlot', 'Add', 'TickLabelInterpreter', 'Latex', 'FontSize', 16, ...
    'Box', 'On', 'LineWidth', 1, 'TickDir', 'both', 'TickLength', [0.005, 0.01]);

h_lines(2) = plot(h_a(1), z_vec*10, E_z{2}, 'LineWidth', 1.5);
h_lines(1) = plot(h_a(1), z_vec*10, ones(size(z_vec)), '-', 'LineWidth', 1.5);
h_gray = plot(h_a(1), z_vec*10, zeros(size(z_vec)), '--', 'LineWidth', 1, 'Color', [1,1,1]*0.5);
uistack(h_gray, 'bottom');

ylabel(h_a(1), 'Normalized co-axial E-field $E_{x}$', 'FontSize', 20, 'Interpreter', 'latex');
xlabel(h_a(1), 'Position along axon $x$ (mm)', 'FontSize', 20, 'Interpreter', 'latex');
xlim(h_a(1), [0, 1]);

h_l = legend(h_a(1), h_lines, {'$E_{\mathrm{uniform}}$', '$E_{\mathrm{micro}}$'});
set(h_l, 'Location', 'Northwest', 'FontSize', 20, 'Box', 'On', 'NumColumns', 1,...
            'Interpreter', 'Latex', 'Interruptible', 'on', 'EdgeColor', [1,1,1], 'Orientation', 'Vertical');

h_lines2(1) = plot(h_a(2), SEinterp_vec, pSE_actual_mean, 'k-', 'LineWidth', 1.5);
patch(h_a(2),[SEinterp_vec; flip(SEinterp_vec)], [pSE_actual_mean-pSE_actual_std; flip(pSE_actual_mean+pSE_actual_std)], ...
             [1,1,1] * 0.5, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

mean_weise = trapz(SEinterp_vec, pSEinterp_vec .* SEinterp_vec);
h_lines2(2) = plot(h_a(2),SEinterp_vec, pSEinterp_vec, 'r-.', 'LineWidth', 1.5);
uistack(h_lines2(2), 'bottom');

ylabel(h_a(2), 'Probablity density function $p(s_{e})$', 'FontSize', 18, 'Interpreter', 'latex');
xlabel(h_a(2), 'E-field ratio $s_{e}=E_{\mathrm{micro}}/E_{\mathrm{uniform}}$', 'FontSize', 18, 'Interpreter', 'latex');
ylim(h_a(2), [0, max(ylim(h_a(2)))]);

plot(h_a(2), [1,1] *mean_weise, ylim(h_a(2)), 'r:', 'LineWidth', 1.5);
plot(h_a(2), [1,1] * mean([mean_x{:}]), ylim(h_a(2)), 'k:', 'LineWidth', 1.5);


h_l = legend(h_a(2), h_lines2, {sprintf('%d E-field profiles, $\\bar{s}_{e}=%.3f$', N,  mean([mean_x{:}])), ...
    sprintf('Weise et al. 2025, $\\bar{s}_{e}=%.3f$',  mean_weise)});
set(h_l, 'Location', 'Northeast', 'FontSize', 20, 'Box', 'On', 'NumColumns', 1,...
            'Interpreter', 'Latex', 'Interruptible', 'on', 'EdgeColor', [1,1,1], 'Orientation', 'Vertical');

save(fullfile('..', 'E-field_profiles.mat'), 'N', 'lambda', 'phi', 'A', 'dz', 'z_vec', 'E_z');
% cm; E-field: normalized -> 1 mV/cm
