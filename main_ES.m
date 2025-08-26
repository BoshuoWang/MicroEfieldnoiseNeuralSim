function results = main_ES_HH(	mod_prmtr, out_ctrl )
%
% mod_prmtr: structure specifying model parameter:
% model_name                Model: 'UF', 'PE'
% model_axon                Axon:  'HH', 'MRG'
% model_size                Axon size: '3um', '0.3um'
% id                        Parameter ID of test case (integer)
%
% out_ctrl: structure specifying output control:
% if_save_data   	Whether to save results in a .mat file for each simulation (logical or 0/1)
% if_write_log  	Whether to write threshold finding process in a .txt log (logical or 0/1)
% if_plot         	Whether to plot threshold finding process (logical or 0/1)

%% Logistics
addpath('Shared functions and data');        % Path of shared functions and data

folder_name = [ mod_prmtr.model_name, '_', mod_prmtr.model_axon, '_', mod_prmtr.model_size];  % Folder for model, add membrane type to folder name          
create_folders(out_ctrl, folder_name);                       % Create subfolders for data, logs, and figures

specify_func = str2func(sprintf('specify_model_ES_%s', mod_prmtr.model_axon));

out_ctrl.log_fid = 0;                                       % Defaults is to display in MATLAB command window
if out_ctrl.if_write_log                                    % Write in .txt file
    logfilename = fullfile(folder_name,'Logs',['log_',num2str(mod_prmtr.id),'.txt']);	% Log filename
    write_fun( out_ctrl.log_fid,{' ', ['Simulation process will be written in: ',logfilename]});
    out_ctrl.log_fid = fopen(logfilename,'w');              % Open log file
end

%% Set parameters for simulation loop
num_sim = 2;                                                % Number of simulations 2
MCE = [0, 0];                                               % 0 using conventional cable equation (CE), 1 using modified CE
field_str = {'macro', 'micro'};                              % 

results = struct(   'th_macro',  NaN, ...                      % Default: straigth axon, CE
                    'th_micro', NaN, 'th_per_diff_micro', NaN...
                );  

%%  Main loop for simulations
t_main = tic;

for ii = 1 : num_sim

    for jj = 1 : 10
        rng("shuffle");
        pause(rand(1));
    end
    
    t_th = tic;
    mod_prmtr.ismicrofield = ii;    % 1: macro; 2: micro
    [solver, stimulation, cable] = specify_func( mod_prmtr );	% Specify parameters for solver, stimulation, and cable

    solver.is_MCE =  MCE(ii);                                           % Whether to use MCE
    
    if ii == 1
        write_fun( out_ctrl.log_fid, solver.txt.log_txt);              	% Output model specification related text
    elseif ~isnan(results.(['th_',field_str{1}]))                          % If threshold found for default simulation (ii == 1, straigth axon, CE), search parameters adjusted
        solver.thresh_find.amp_init = results.(['th_',field_str{1}]) * (0.9 + randn(1) * 0.01);	% Initial search amplitude set to 80% of default threshold
        solver.thresh_find.factor   = (solver.thresh_find.factor)^(1/2);                        % Smaller factor more robust search
    end
        
    write_fun(out_ctrl.log_fid, {'-----------------------------------------------------------',...
                                 ['Solver: ',field_str{ii}],' '});

    if out_ctrl.if_plot         % Set up figure
        out_ctrl.h_fig = figure('Position',[00 00 1400 800],'Color',[1,1,1]);
        axes('position',[0.05,0.95,0.9,0.00]);box off; axis off;
        title(  [solver.txt.fig_title, field_str{ii}], 'Interpreter','latex','FontSize',14);
    end
    
    results.(['th_',field_str{ii}]) = threshold_finding( solver, stimulation, cable, out_ctrl);

    if ii ~= 1                                                  % Percentage difference of threshold compared to default
        if ~isnan(results.(['th_',field_str{1}]))
            results.(['th_per_diff_',field_str{ii}]) = ( results.(['th_',field_str{ii}]) / results.(['th_',field_str{1}]) - 1) * 100;
            
        elseif ~isnan(results.(['th_',field_str{ii}]))           % If no threshold obtained for default
            results.(['th_per_diff_',field_str{ii}]) = -100;     % Then any obtained threshold is 100% less
        end
        write_fun(out_ctrl.log_fid,	{sprintf('Percentage difference: %2.3f %%.', results.(['th_per_diff_',field_str{ii}]))});    % Write/display results
    end
                      
    if out_ctrl.if_plot                                         % Save and close figures
        figure_filename = fullfile(folder_name,'Figures', ['Fig_',num2str(mod_prmtr.id),'_',field_str{ii},'.fig']);
        % saveas(out_ctrl.h_fig,figure_filename, 'fig');
        im = frame2im(getframe(out_ctrl.h_fig));
        imwrite( im, [figure_filename(1:end-4), '.tif'], 'tif', 'WriteMode', 'overwrite', 'Resolution', 300, 'Compression', 'lzw');

        close(out_ctrl.h_fig);
    end
    
    write_fun(out_ctrl.log_fid, {' ', sprintf('Run time for search: %3.2f min.',toc(t_th)/60), ' '});
end

%% Saving results and closing files
write_fun( out_ctrl.log_fid,    {'-----------------------------------------------------------',...
                                 sprintf('Total run time: %3.3f min.', toc(t_main)/60)});    % Output run time 
if out_ctrl.if_save_data
    filename = fullfile(folder_name,'Results',['result_',num2str(mod_prmtr.id),'.mat']);
    save(filename,'mod_prmtr', 'results');
    write_fun(out_ctrl.log_fid, {sprintf('Results saved in %s.', filename)});
end
if out_ctrl.if_write_log
    out_ctrl.log_fid = fclose(out_ctrl.log_fid);
end
rmpath('Shared functions and data');
end
