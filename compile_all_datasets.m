
model = {'UF_HH_3um', 'UF_HH_0.3um', 'PE_HH_3um', 'PE_HH_0.3um', ...
         'UF_RMG_3um', 'UF_RMG_0.3um', 'PE_RMG_3um', 'PE_RMG_0.3um'};

for ii = 1 : length(model)
    model_name = model{ii};
    compile(model_name);
end


function compile(model_name)
PW_vec = [ kron([1e-2,1e-1,1e0], [1,3.1]), 1e1];
% PW in ms, 1 us to 10 ms; >10 ms rheobase (steady state for linear
% membrane); 25 points (4 decades)

if strncmp(model_name,'UF',2)                  % Uniform Field
    N_vec = 1:100;
    [NN, PPWW] = ndgrid(N_vec, PW_vec);
else                                            % Electrode stimulation
    N_vec = 1:40;
    [NN, PPWW] = ndgrid(N_vec, PW_vec);
end

NaN_matrix = NaN(size(PPWW));

compiled_results = struct(  'model',model_name,...
    'N',NN,'PW',PPWW,...
    'parameter_id',NaN_matrix...
    );


for ii = 1: numel(PPWW)
    filename = fullfile(model_name,'Results',['result_',num2str(ii),'.mat']);
    fprintf('Loading file %s\n', filename)
    if exist(filename,'file') > 0
        load(filename,'mod_prmtr', 'results');

        parameter_id = mod_prmtr.id;
        compiled_results.parameter_id(parameter_id)	= parameter_id;
        compiled_results.th_macro(parameter_id) = results.th_macro;
        compiled_results.th_micro(parameter_id) = results.th_micro;
        compiled_results.th_per_diff_micro(parameter_id) = results.th_per_diff_micro;

    else
        disp(['File does not exist: ',filename]);
    end
end

filename = fullfile(model_name,[model_name,'_compiled_result.mat']);
save(filename,'compiled_results');
disp(['Compiled ', model_name,'.']);
end
