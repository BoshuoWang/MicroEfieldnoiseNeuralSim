function create_folders(out_ctrl, folder_name)

% Create subfolders, if they not already exist
% For saving data
if out_ctrl.if_save_data
    subfoldername = fullfile(folder_name, 'Results');
    if ~exist(subfoldername,'dir')
        mkdir(subfoldername);
    end
end

% For logs
if out_ctrl.if_write_log
    subfoldername = fullfile(folder_name, 'Logs');
    if ~exist(subfoldername,'dir')
        mkdir(subfoldername);
    end
end

% For figures
if out_ctrl.if_plot
    subfoldername = fullfile(folder_name, 'Figures');
    if ~exist(subfoldername, 'dir')
        mkdir(subfoldername);
    end
end

end