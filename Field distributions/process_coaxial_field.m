% clearvars;
opt = optimset('fminsearch');                   % parameters for curve fitting function
opt = optimset(opt,'MaxFunEvals',10000,'MaxIter',10000,'TolX',1e-10,'TolFun',1e-8,'Display','off');
dst = 5;

AF_method = 2; % 1: COMSOL, 2: dEx/dx, 3: d^2V/dx^2

lambda_axon = 15e-6;
fc = 1/(2*pi*lambda_axon);
fs = (1e-6)^-1;

[b, a] = butter(2, fc/(fs/2));

lambda_axon2 = 15e-6;
fc2 = 1/(2*pi*lambda_axon2);
fs2 = fs/dst;

[b2, a2] = butter(2, fc2/(fs2/2));

%%
filename = 'homogeneous_y118z81';
data1 = importfile_coax([filename, '.txt']);

ind = ~isnan(data1(: ,7));

x1 = data1(ind,1);
Ex1 = data1(ind,4);
Ey1 = data1(ind,5);
Ez1 = data1(ind,6);
V1 = data1(ind,7);
AF1 = data1(ind,8);
V1_max = max(V1);

x1_AF = x1(1:dst:end);
x1_AF =[-flipud(x1_AF); x1_AF(2:end)];

x1 = [-flipud(x1); x1(2:end)];
V1 = [flipud(V1); V1(2:end)];
Ex1 = [-flipud(Ex1); Ex1(2:end)];
Ey1 = [flipud(Ey1); Ey1(2:end)];
Ez1 = [flipud(Ez1); Ez1(2:end)];
AF1 = [flipud(AF1); AF1(2:end)];

filename = 'micro3D_y118z81';
data2 = importfile_coax([filename, '.txt']);

ind = ~isnan(data2(: ,7));

x2 = data2(ind,1);
Ex2 = data2(ind,4);
Ey2 = data2(ind,5);
Ez2 = data2(ind,6);
V2 = data2(ind,7);
AF2 = data2(ind,8);
V2_max = max(V2);

x2_AF = x2(1:dst:end);
x2_AF =[-flipud(x2_AF); x2_AF(2:end)];

x2 = [-flipud(x2); x2(2:end)];
V2 = [flipud(V2); V2(2:end)];
Ex2 = [-flipud(Ex2); Ex2(2:end)];
Ey2 = [flipud(Ey2); Ey2(2:end)];
Ez2 = [flipud(Ez2); Ez2(2:end)];
AF2 = [flipud(AF2); AF2(2:end)];


Ex2_filt = filtfilt(b, a, Ex2);

[~, ind_c] = min(abs(x1));
ind = [flip(ind_c:-dst:1), ind_c+dst:dst:length(x2)];
        
switch AF_method 
    case 1
        AF1 = AF1(ind);
        AF2 = AF2(ind);
        % AF2_filt = filtfilt(b2, a2, AF2);
    case 2
        tmp = Ex1(ind);
        AF1 = ([0;diff(tmp)]+[diff(tmp);0])/2/(x1(dst+1)-x1(1));
        tmp = Ex2(ind);
        AF2 = ([0;diff(tmp)]+[diff(tmp);0])/2/(x2(dst+1)-x2(1));

        tmp = Ex2_filt(ind);
        % AF2_filt = ([0;diff(tmp)]+[diff(tmp);0])/2/(x2(dst+1)-x2(1));
        
    case 3
        AF1 = -diff(V1(ind),2)/(x1(dst+1)-x1(1))^2;
        AF2 = -diff(V2(ind),2)/(x2(dst+1)-x2(1))^2;
        % AF2_filt = filtfilt(b2, a2, AF2);
end

% AF1_filt = filtfilt(d, AF1);
AF2_filt = filtfilt(b2, a2, AF2);
% AF2_filt = movmean(AF2, 3);


lambdaVmax = V1_max/V2_max;
[lambdaE, ~, exit_flag] = fminsearch(@(x) JError(x, Ex2, Ex1), 1, opt);
[lambda_fitE2, fit, exit_flag] = fminsearch(@(x) JError(x, [Ex2(:), Ey2(:), Ez2(:)], [Ex1(:), Ey1(:), Ez1(:)], 1), 1, opt);

[lambdaV, ~, exit_flag] = fminsearch(@(x) JError(x, V2, V1), 1, opt);

lambda = [lambdaVmax, lambdaV, lambdaE, lambda_fitE2]
max(AF2)/max(AF1)
%%
clearvars h*
figure('Color', 'w', 'Position', [20, 20, 2100, 600]);
col = lines(7);

h_a(1) = subplot(1,3,1);
h_a(2) = subplot(1,3,2);
h_a(3) = subplot(1,3,3);

set(h_a, 'NextPlot', 'Add', 'TickLabelInterpreter', 'Latex', 'FontSize', 14, ...
    'Box', 'Off', 'LineWidth', 1, 'TickDir', 'both', 'TickLength', [0.005, 0.01], ...
    'XMinorTick', 'off', 'XMinorGrid', 'off', ...
    'YGrid', 'off', 'YMinorTick', 'off', 'YMinorGrid', 'off');

xlabel(h_a, '$x$ (mm)', 'FontSize', 20, 'Interpreter', 'latex');
ylabel(h_a(1), '$V \ \mathrm{(mV)}$ ', 'FontSize', 20, 'Interpreter', 'latex');
ylabel(h_a(2), '$E_x \ \mathrm{(mV/mm)}$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel(h_a(3), '$\mathrm{d^2}V / {\mathrm{d}x}^2 \ \mathrm{(mV/{mm}^2)}$', 'FontSize', 20, 'Interpreter', 'latex');

set(h_a, 'XLim', [-0.5, 0.5]);


col_ind_2 = 6;
col_ind_filt = 1;
col_ind_1 = 2;

plot(h_a(1), x2/1e-3, -V2/1e-3, 'Color', col(col_ind_2,:), 'LineWidth', 2.5);
plot(h_a(2), x2/1e-3, -Ex2, 'Color', col(col_ind_2,:),  'LineWidth', 1);


% f = lambdaV;
f = 0.59;

plot(h_a(1), x1/1e-3, -V1/1e-3/f, 'Color', col(col_ind_1,:), 'LineWidth', 1.5);
plot(h_a(2), x1/1e-3, -Ex1/f, 'Color', col(col_ind_1,:), 'LineWidth', 1.5);

plot(h_a(1), x1/1e-3, -V1/1e-3, 'k--', 'LineWidth', 1.5);
plot(h_a(2), x1/1e-3, -Ex1, 'k--', 'LineWidth', 1.5);

plot(h_a(2), x2/1e-3, -Ex2_filt, '--', 'Color', col(col_ind_filt,:), 'LineWidth', 1.5);


switch AF_method 
    case 1
        plot(h_a(3), x2_AF/1e-3, AF2*1e-3, 'Color', col(col_ind_2,:),  'LineWidth', 1);
        plot(h_a(3), x1_AF/1e-3, AF1/f*1e-3, 'Color', col(col_ind_1,:),  'LineWidth', 1.5);
        plot(h_a(3), x1_AF/1e-3, AF1*1e-3, 'k--', 'LineWidth', 1.5);
        h_filt = plot(h_a(3), x2_AF/1e-3, AF2_filt*1e-3, '--', 'Color', col(col_ind_filt,:), 'LineWidth', 1.5);        
        title_str = 'COMSOL';
    case 2
        plot(h_a(3), x2_AF/1e-3, AF2*1e-3, 'Color', col(col_ind_2,:),  'LineWidth', 1);
        plot(h_a(3), x1_AF/1e-3, AF1/f*1e-3, 'Color', col(col_ind_1,:),  'LineWidth', 1.5);
        plot(h_a(3), x1_AF/1e-3, AF1*1e-3, 'k--', 'LineWidth', 1.5);
        h_filt = plot(h_a(3), x2_AF/1e-3, AF2_filt*1e-3, '--', 'Color', col(col_ind_filt,:), 'LineWidth', 1.5);
        
        title_str = ['MATLAB $\mathrm{d}E_{\mathrm{x}} / {\mathrm{d}x}$: ', num2str((x1(dst+1)-x1(1))/1e-6), ' $\mu m$'];
    case 3
        plot(h_a(3), x2_AF(2:end-1)/1e-3, AF2*1e-3, 'Color', col(col_ind_2,:), 'LineWidth', 1);
        plot(h_a(3), x1_AF(2:end-1)/1e-3, AF1/f*1e-3, 'Color', col(col_ind_1,:), 'LineWidth', 1.5);
        plot(h_a(3), x1_AF(2:end-1)/1e-3, AF1*1e-3, 'k--', 'LineWidth', 1.5);
        h_filt = plot(h_a(3), x2_AF(2:end-1)/1e-3, AF2_filt*1e-3, '--', 'Color', col(col_ind_filt,:), 'LineWidth', 1.5);
        
        title_str = ['MATLAB $\mathrm{d^2}V / {\mathrm{d} x}^2$: ', num2str((x1(dst+1)-x1(1))/1e-6), ' $\mu m$'];
end

set(h_a(1), 'YLim', h_a(1).YLim * 1.25);

h_leg(1) = legend(h_a(1), ...
    {'Microscopic ($\sigma$ = 1 S/m)', ...
     ['Homogenous, corrected ($\sigma$ = ', sprintf('%1.1f', f) ,' S/m)'], ...
      'Homogenous, uncorrected ($\sigma$ = 1 S/m)'});


h_leg(2) = legend(h_a(3), h_filt, ...
    {'Microscopic, filtered ($\sigma$ = 1 S/m)', ...
     });

set(h_leg, 'Location', 'Southwest', 'FontSize', 15, 'Box', 'Off', ...
    'Interpreter', 'Latex', 'Interruptible', 'on', 'EdgeColor', [1,1,1]*0.5,...
    'Orientation', 'Vertical', 'NumColumns', 1);


%%
load("pSE_distribution.mat");

pSE_pointsource = Ex2./(Ex1/f);
SE_vec = (0:0.1:5)';
SEint_vec = (SE_vec(1:end-1) + SE_vec(2:end))/2;

ind_end = length(x1);
ind_mid = (ind_end+1)/2;
ind_vec = ind_mid : 23 : ind_end;

N = length(ind_vec);
for nn = N : -1 : 1
    pSE_actual{nn} = histcounts(abs(pSE_pointsource(1 : ind_vec(nn))), SE_vec, 'Normalization', 'pdf')';
    
    mean_x{nn} = trapz(SEint_vec, SEint_vec .* pSE_actual{nn});
    std_x{nn} = sqrt( trapz(SEint_vec, (SEint_vec - mean_x{nn}).^2 .* pSE_actual{nn}) );
    sk_x{nn} = trapz(SEint_vec, ( (SEint_vec - mean_x{nn}) / std_x{nn} ).^3 .* pSE_actual{nn});
end
pSE_actual_mean = mean([pSE_actual{:}], 2);
pSE_actual_std = std([pSE_actual{:}], 0, 2);

col = lines(7);
clearvars h*
for hh = 1 %: 2
    h_f(hh) = figure('Color', 'w', 'Position', [75 * hh, 50, 640, 360]);
    h_a(hh) = axes(h_f(hh));
    h_ax_1(hh) = axes(h_f(hh), 'Position',[0.005, 0.945, 0.05, 0.05],'Clipping','off','Box','off');
    axis(h_ax_1(hh), 'off');
    axis(h_ax_1(hh), [0, 1, -1, 0]);
end
text(h_ax_1(1), 0, 0, 0, '$\mathbf{B_{ii}}$', 'Interpreter', 'Latex', 'FontSize', 30, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
h_f(2) = openfig(fullfile('..', 'Figures_Full', '2C_tmp.fig'));
h_a(2) = h_f(2).Children(end);
h_f(3) = openfig(fullfile('..', 'Figures_Full', '2D_tmp.fig'));
h_a(3) = h_f(3).Children(end);

set(h_a, 'NextPlot', 'Add', 'TickLabelInterpreter', 'Latex', 'FontSize', 16, ...
    'Box', 'On', 'LineWidth', 1, 'TickDir', 'both', 'TickLength', [0.005, 0.01]);

xlabel(h_a(1), 'Position along axon $x$ (mm)', 'FontSize', 18, 'Interpreter', 'latex');
ylabel(h_a(1), '$E_x \ \mathrm{(mV/mm)}$', 'FontSize', 18, 'Interpreter', 'latex');

set(h_a(1), 'XLim', [-0.5, 0.5]);
h_lines(2) = plot(h_a(1), x2/1e-3, -Ex2, 'LineWidth', 1.5, 'Color', col(2, :));
h_lines(1) = plot(h_a(1), x1/1e-3, -Ex1/f, 'LineWidth', 1.5, 'Color', 'k');
h_gray = plot(h_a(1), x1/1e-3, zeros(size(x1)), '--', 'LineWidth', 1, 'Color', [1,1,1]*0.5);
uistack(h_gray, 'bottom');

h_l = legend(h_a(1), h_lines, {'$E_{\mathrm{macro}}$', '$E_{\mathrm{total}}$'});
set(h_l, 'Location', 'Northeast', 'FontSize', 18, 'Box', 'On', 'NumColumns', 1,...
    'Interpreter', 'Latex', 'Interruptible', 'on', 'EdgeColor', [1,1,1], 'Orientation', 'Vertical');
ylim(h_a(1), [-30,30]);


h_lines2(1) = plot(h_a(2), SEint_vec, pSE_actual_mean, '-', 'LineWidth', 1.5, 'Color', col(2, :));
patch(h_a(2), [SEint_vec; flip(SEint_vec)], [pSE_actual_mean-pSE_actual_std; flip(pSE_actual_mean+pSE_actual_std)], ...
             col(2, :), 'FaceAlpha', 0.5, 'EdgeColor', 'none');


mean_weise = trapz(SEinterp_vec, pSEinterp_vec .* SEinterp_vec);

ylabel(h_a(2), '$p(s_{e})$', 'FontSize', 18, 'Interpreter', 'latex');
xlabel(h_a(2), 'E-field ratio $s_{e}=E_{\mathrm{total}}/E_{\mathrm{macro}}$', 'FontSize', 18, 'Interpreter', 'latex');
set(h_a(2), 'YLim', [0,2])
ylim(h_a(2), [0, max(ylim(h_a(2)))]);

h_a(2).Children = h_a(2).Children([1,3,2,4,5]);


plot(h_a(2), [1,1] * mean_weise, ylim(h_a(2)), 'r:', 'LineWidth', 1.5);
plot(h_a(2), [1,1] * mean([mean_x{:}]), ylim(h_a(2)), ':', 'LineWidth', 1.5, 'Color', 'k');


h_l = legend(h_a(2), h_a(2).Children([6,5,7]), ...
    {h_f(2).Children(2).String{1}, ...
    sprintf('Spherical electrode, $\\bar{s}_{e}=%.3f$', mean([mean_x{:}])), ...
    h_f(2).Children(2).String{2}});
set(h_l, 'Location', 'Northeast', 'FontSize', 16, 'Box', 'On', 'NumColumns', 1,...
            'Interpreter', 'Latex', 'Interruptible', 'on', 'EdgeColor', [1,1,1], 'Orientation', 'Vertical');

dz = median(diff(x1)) / 1e-3;
dz = dz/5;
E_noise = Ex2-Ex1/f;
E_noise = resample(E_noise, 5, 1);

fs = 1/(dz);                                      % mm^-1, sampling frequency
NFFT = 2^nextpow2(length(E_noise));
df = fs/NFFT;                                   % mm^-1, frequency resolution
freq_vec = 0 : df : fs- df;               % frequency vector of FFT, only cover DC to fs/2
ind_freq = 1 : (NFFT/2+1);

FFT_tmp = fft(E_noise, NFFT)/sqrt(NFFT);               % NFFT scaling is to adjust for proper amplitude
FFT_vec = FFT_tmp(ind_freq);                             % discard mirrowed
psd_vec = abs(FFT_vec).^2;

FFT_vec = FFT_vec/max(abs(FFT_vec));

plot(h_a(3), freq_vec(ind_freq), abs(FFT_vec), 'Color', col(2, :));
ylabel(h_a(3), {'Normalized spectrum'}, 'FontSize', 18, 'Interpreter', 'latex');
xlabel(h_a(3), 'Spatial frequency $(\mathrm{mm^{-1}})$', 'FontSize', 18, 'Interpreter', 'latex');
set(h_a(3), 'XScale', 'log', 'YScale', 'log', 'xlim', [1, 1000], 'YLim', [0.99e-5, 1e0], 'YTick', 10.^(-6:0));


h_a(1).Position([1,3]) = [0.15, 0.8];
h_a(2).Position([1,3]) = [0.15, 0.8];
h_a(3).Position([1,3]) = [0.15, 0.8];

saveas(h_f(1), fullfile('..', 'Figures_Full', '2Bii.fig'));
saveas(h_f(2), fullfile('..', 'Figures_Full', '2C.fig'));
saveas(h_f(3), fullfile('..', 'Figures_Full', '2D.fig'));

%%
save(fullfile('..', 'PointSource.mat'), 'x1', 'x2', 'V1', 'V2', 'Ex1', 'Ex2', 'f');
% m, V, V/m, for 8 μA