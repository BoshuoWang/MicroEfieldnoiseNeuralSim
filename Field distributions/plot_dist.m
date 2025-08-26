clearvars
load(fullfile('..', 'E-field_profiles.mat'), 'N', 'lambda', 'phi', 'A', 'dz', 'z_vec', 'E_z');
load("pSE_distribution.mat");
dSE = mean(diff(SE_vec));
SE_vec2 = [(-1 : dSE : -dSE)' ; SE_vec];
SEinterp_vec2 = (SE_vec2(1:end-1) + SE_vec2(2:end))/2;


for nn = N : -1 : 1
    pSE_actual{nn} = histcounts(abs(E_z{nn}), SE_vec, 'Normalization', 'pdf')';
    pSE_actual_signed{nn} = histcounts(E_z{nn}, SE_vec2, 'Normalization', 'pdf')';
    
    mean_x{nn} = trapz(SEinterp_vec, SEinterp_vec .* pSE_actual{nn});
    std_x{nn} = sqrt( trapz(SEinterp_vec, (SEinterp_vec - mean_x{nn}).^2 .* pSE_actual{nn}) );
    sk_x{nn} = trapz(SEinterp_vec, ( (SEinterp_vec - mean_x{nn}) / std_x{nn} ).^3 .* pSE_actual{nn});
end


pSE_actual_mean = mean([pSE_actual{:}], 2);
pSE_actual_std = std([pSE_actual{:}], 0, 2);

pSE_actual_signed_mean = mean([pSE_actual_signed{:}], 2);
pSE_actual_signed_std = std([pSE_actual_signed{:}], 0, 2);


%%
col = lines(7);
   
clearvars h*
for hh = 1 : 3
    h_f(hh) = figure('Color', 'w', 'Position', [75 * hh, 50, 640, 360]);
    h_a(hh) = axes(h_f(hh));
    h_ax_1(hh) = axes(h_f(hh), 'Position',[0.005, 0.945, 0.05, 0.05],'Clipping','off','Box','off');
    axis(h_ax_1(hh), 'off');
    axis(h_ax_1(hh), [0, 1, -1, 0]);
end

set(h_a, 'NextPlot', 'Add', 'TickLabelInterpreter', 'Latex', 'FontSize', 16, ...
    'Box', 'On', 'LineWidth', 1, 'TickDir', 'both', 'TickLength', [0.005, 0.01]);

h_lines(2) = plot(h_a(1), z_vec*10, E_z{2}, 'LineWidth', 1.5, 'Color', col(1,:));
h_lines(1) = plot(h_a(1), z_vec*10, ones(size(z_vec)), '-', 'LineWidth', 1.5, 'Color', 'k');
h_gray = plot(h_a(1), z_vec*10, zeros(size(z_vec)), '--', 'LineWidth', 1, 'Color', [1,1,1]*0.5);
uistack(h_gray, 'bottom');

ylabel(h_a(1), {'Normalized $E_{x}$'}, 'FontSize', 18, 'Interpreter', 'latex');
xlabel(h_a(1), 'Position along axon $x$ (mm)', 'FontSize', 18, 'Interpreter', 'latex');
xlim(h_a(1), [0, 1]);

h_l = legend(h_a(1), h_lines, {'$E_{\mathrm{macro}}$', '$E_{\mathrm{total}}$'});
set(h_l, 'Location', 'Northwest', 'FontSize', 18, 'Box', 'On', 'NumColumns', 1,...
            'Interpreter', 'Latex', 'Interruptible', 'on', 'EdgeColor', [1,1,1], 'Orientation', 'Vertical');
ylim(h_a(1), [-0.7, 2.8]);

h_lines2(1) = plot(h_a(2), SEinterp_vec, pSE_actual_mean, '-', 'LineWidth', 1.5, 'Color', col(1, :));
patch(h_a(2),[SEinterp_vec; flip(SEinterp_vec)], [pSE_actual_mean-pSE_actual_std; flip(pSE_actual_mean+pSE_actual_std)], ...
             col(1, :), 'FaceAlpha', 0.35, 'EdgeColor', 'none');

mean_weise = trapz(SEinterp_vec, SEinterp_vec .* pSEinterp_vec);
std_weise = sqrt(trapz(SEinterp_vec, (SEinterp_vec - mean_weise).^2 .* pSEinterp_vec));
sk_weise = trapz(SEinterp_vec, ((SEinterp_vec - mean_weise)/std_weise).^3 .* pSEinterp_vec);

h_lines2(2) = plot(h_a(2),SEinterp_vec, pSEinterp_vec, 'r-.', 'LineWidth', 1.5);
uistack(h_lines2(2), 'bottom');

ylabel(h_a(2), '$p(s_{e})$', 'FontSize', 18, 'Interpreter', 'latex');
xlabel(h_a(2), 'E-field ratio $s_{e}=E_{\mathrm{total}}/E_{\mathrm{macro}}$', 'FontSize', 18, 'Interpreter', 'latex');
ylim(h_a(2), [0, 1.1]);


h_l = legend(h_a(2), h_lines2, {sprintf('%d E-field profiles, $\\bar{s}_{e}=%.3f$', N,  mean([mean_x{:}])), ...
    sprintf('Weise et al. 2025, $\\bar{s}_{e}=%.3f$',  mean_weise)});
set(h_l, 'Location', 'Northeast', 'FontSize', 16, 'Box', 'On', 'NumColumns', 1,...
            'Interpreter', 'Latex', 'Interruptible', 'on', 'EdgeColor', [1,1,1], 'Orientation', 'Vertical');


fs = 1/(dz*10);                                      % mm^-1, sampling frequency
NFFT = 2^nextpow2(length(E_z{2}));
df = fs/NFFT;                                   % mm^-1, frequency resolution
freq_vec = 0 : df : fs- df;               % frequency vector of FFT, only cover DC to fs/2
ind_freq = 1 : (NFFT/2+1);

FFT_tmp = fft(E_z{2}-mean(E_z{2}), NFFT)/sqrt(NFFT);               % NFFT scaling is to adjust for proper amplitude
FFT_vec = FFT_tmp(ind_freq);                             % discard mirrowed
psd_vec = abs(FFT_vec).^2;

FFT_vec = FFT_vec/max(abs(FFT_vec));

plot(h_a(3), freq_vec(ind_freq), abs(FFT_vec), 'Color', col(1, :));
ylabel(h_a(3), 'Spectral density $(\mathrm{a.u.})$', 'FontSize', 18, 'Interpreter', 'latex');
xlabel(h_a(3), 'Spatial frequency $(\mathrm{mm^{-1}})$', 'FontSize', 18, 'Interpreter', 'latex');
set(h_a(3), 'XScale', 'log', 'YScale', 'log', 'xlim', [1, 2000]);

h_a(1).Position([1,3]) = [0.15, 0.8];
h_a(2).Position([1,3]) = [0.15, 0.8];
h_a(3).Position([1,3]) = [0.15, 0.8];

%%
text(h_ax_1(1), 0, 0, 0, '$\mathbf{A_{ii}}$', 'Interpreter', 'Latex', 'FontSize', 30, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
text(h_ax_1(2), 0, 0, 0, '$\mathbf{C}$', 'Interpreter', 'Latex', 'FontSize', 30, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
text(h_ax_1(3), 0, 0, 0, '$\mathbf{D}$', 'Interpreter', 'Latex', 'FontSize', 30, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');

if ~exist(fullfile('..', 'Figures_Full'), 'dir')
    mkdir(fullfile('..', 'Figures_Full'));
end

saveas(h_f(1), fullfile('..', 'Figures_Full', '2Aii.fig'));
saveas(h_f(2), fullfile('..', 'Figures_Full', '2C_tmp.fig'));
saveas(h_f(3), fullfile('..', 'Figures_Full', '2D_tmp.fig'));