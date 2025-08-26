function [solver, stimulation, cable] = specify_model_ES_HH(modprmtr)

PW_vec = [ kron([1e-2, 1e-1, 1e0], [1, 3.1]), 1e1];
% PW in ms, 10 us to 10 ms; >10 ms rheobase (steady state for linear
% membrane)

if strncmp(modprmtr.model_name, 'UF', 2)                  % Uniform Field
    mod_str = 'Axon terminal';
    load('E-field_profiles.mat', 'N', 'lambda', 'phi', 'A', 'dz', 'z_vec', 'E_z');
    % cm; E-field: normalized -> 1 mV/cm
    
    dz_vec = dz;
    
    L_straight = max(z_vec);         % Axon length, in cm
    modprmtr.mod = 1;
    
    N_vec = (1:100)';
    [NN, PPWW] = ndgrid(N_vec, PW_vec);
    
    N = NN(modprmtr.id);
    PW = PPWW(modprmtr.id);
    log_str = sprintf('N:\t%d', N);
    title_str = ['$$ N=', num2str(N,'%d'),'$$; '];
    dt = 1e-3;              % Time step 2 us, in ms
else                                            % electrode stimulation
    modprmtr.isdisk = 0;
    mod_str = 'Point electrode';
    
    load('PointSource.mat', 'x1', 'x2', 'V1', 'V2', 'Ex1', 'Ex2', 'f');
    V1 = V1 * 1e3;  % V -> mV
    V2 = V2 * 1e3;
    x1 = x1 * 1e2;  % m -> cm
    x2 = x2 * 1e2;

    dz = median(diff(x1));

    modprmtr.mod = 3;
    
    N_vec = 1:40;
    [NN, PPWW] = ndgrid(N_vec, PW_vec);
    
    N = NN(modprmtr.id);
    PW = PPWW(modprmtr.id);
    log_str = sprintf('N:\t%d', N);
    title_str = ['$$ N=', num2str(N,'%d'),'$$; '];
    dt = 1e-3;              % Time step 2 us, in ms
end

solver.txt.log_txt = {  sprintf('Model type: \t%s', mod_str),...
                        sprintf('Parameter ID:\t%d',modprmtr.id), ...
                        sprintf('Pulse duration in ms:\t%2.4f',PW),...
                        log_str};
solver.txt.fig_title = [mod_str,'; ',...
                        'Parameter $$ID=',num2str(modprmtr.id),'$$; ',...
                        '$$ PW = ',num2str(PW,'%3.3f'),'\: \rm{ms}$$; ',...
                        title_str];

%% Cellular parameters
% Specific to C & PS model (Roth & Basser 1990; Schnabel & Johannes, 2001;
% Neu, 2016a,b)
T = 23.5;       % degree Celcius; Roth & Basser 1990: 18.5 C

% Conductivities of extra- and intra-cellular spaces
sigma_i = 28.2 ;            % Intracellular conductivity, in mS/cm;
sigma_e = 10;               % Extracellular conductivity, in mS/cm;

if strncmp(modprmtr.model_size, '3um', 3)
    R_a = 3e-4;              % Axon radius, in cm; 3 um radius
elseif strncmp(modprmtr.model_size, '0.3um', 5)
    R_a = 0.3e-4;              % Axon radius, in cm; 0.3 um radius
end
c_m = 1;                    % Membrane capacitance, uF/cm^2
% Cellular time constant for cylindrical cell is 13.55 ns per um in radius;
% tau_c = R * c_m * (sigma_i^-1 + sigma_e^-1)
% 3 um radius, tau_c =  40.6 ns; 
% 0.1 us = 2.5*tau_c    ->  reaching 91.8% of IP
% 0.2 us = 5 tau_c      ->  reaching 99.3% of IP

% HH parameters: original Hodgkin-Huxley
V_rest = -70;   % mV; Roth & Basser 1990: -65 mV; Cartee 2000, Rattay & Aberham 1993: -70 mV
E_Na = V_rest + 115;        g_Na = 120;     % mV & mS/cm^2;
E_K  = V_rest - 12;         g_K  = 36;   
E_L  = V_rest + 10.6;       g_L  = 0.3;    
% g_bar = 0.6773 mS/cm^2 at rest -> r_m = 1.476 kOhm*cm^2 


% Calculated parameters
freq = 100e-3;                                          % 100 Hz in kHz
lambda_100 = sqrt(R_a * sigma_i / (2*pi*freq*c_m));     % in cm; lambda(w) = lambda_DC * sqrt( 2 / (1 + sqrt( 1 + (2*w*c_m/(sigma_i*R))^2 * lambda_DC^4 ) ) )
d_lambda = 0.1;                                         % Interval is determined by d_lambda rule
dz = lambda_100 * d_lambda;                             % For R = 3um, lambda100 = 1.16 mm
                                                          % d_lambda = 0.3 -> dz = 348 um; d_lambda = 0.1 -> dz = 116 um
dz = 2e-4; %    2 um in cm
[m,~,h,~,n,~] = HH(V_rest, T, V_rest);                  % Resting ion channel parameters
g_m_rest = g_Na * m^3 * h + g_K * n^4 + g_L;            % g_m is 0.6773 mS/cm^2 at rest
lambda_DC = sqrt(R_a * sigma_i / g_m_rest /2);          % DC length constant, for R = 3 um, lambda_DC = 0.79 mm

% lambda_DC = sqrt(R_a * sigma_i / g_L /2);           
lambda_w = lambda_DC * sqrt( 2 / (1 + sqrt( 1 + (2*(2*pi*freq)*c_m/(sigma_i*R_a))^2 * lambda_DC^4 ) ) );
tau = c_m / g_m_rest;

% Emperical parameter from simulation
v = 0.2335;                   % cm/ms; conduction speed 2.335 mm/ms
v = v * (R_a/3e-4); 

%% Specify cable
N_theta = 15;               % Discretization of azimuthal angle
d_theta =  pi / N_theta;    % Interval for integration, in radian
theta = linspace( d_theta/2, pi - d_theta/2, N_theta );         % Integration points between 0 and pi, row vector

switch modprmtr.mod
    case 1
        axon_length = L_straight;       % Axon length, 1 cm in cm
        N_comp = ceil(axon_length / dz);
    
        axon_length = N_comp * dz;
        z = linspace( dz/2, axon_length-dz/2, N_comp )';   % axial coordinates, in cm. Row vector
        dz = dz * ones(N_comp,1);
    
    case 3
        axon_neg_length = max(x1);
        axon_pos_length = axon_neg_length*(N-1)/40;
        axon_neg_length = axon_neg_length + 0.2;  % +2 mm
        
        N_neg_comp_uni = floor(axon_neg_length/dz);
        N_pos_comp_uni = floor(axon_pos_length/dz);
        
        axon_neg_length =  - N_neg_comp_uni * dz;
        axon_pos_length = N_pos_comp_uni * dz;
        
        axon_length = axon_neg_length + axon_pos_length;
        N_comp = N_neg_comp_uni + N_pos_comp_uni+1;
        z = linspace( axon_neg_length, axon_pos_length, N_comp )';     % axial coordinates, in cm. Row vector
        dz = dz * ones(N_comp,1);
end
          
R_i = dz / (sigma_i*pi*R_a^2);       % Axial resistance between nodes, in kOhm
Area = 2 * pi * R_a * dz;        % Element area; cm^2
C_m = c_m * Area;                   % Node membrane capacitance, in uF
ones_A = ones(N_comp, 1);                     % Column vector; empty array
TP_weight = ones(N_comp, N_theta) * d_theta / pi;
cable = struct( 'N_comp',N_comp,...         % Number of compartments
                'z',z,...                   % Center coordinates of compartment, cm (along local longitudinal axis)
                'dz',dz,...                 % Compartment length, cm
                'R', R_a * ones_A,...       % Compartment radius, cm
                'R_i',R_i,...               % Compartment axial resistance, kOhm
                'R_i_left',R_i,...          % Axial resistance to left neighbor, kOhm
                'R_i_right',R_i,...         % Axial resistance to right neighbor, kOhm
                'Area',Area,...             % Compartment Area, cm^2
                'C_m',C_m,...               % Compartment capacitance, uF
                'TP_dim', ones_A,...        % Dimension: 1 for cylindrical; 2 for sphercial
                'TP_weight',TP_weight...    % Integration weights, row vector for each compartment
                );

% Axial resistance; sealed ends are reflected by d_phi_e_left/right = 0 at terminals
cable.R_i_left(1) = inf; 
cable.R_i_left(2:end)    = ( cable.R_i(1:end-1) + cable.R_i(2:end))/2;
cable.R_i_right(1:end-1) = ( cable.R_i(1:end-1) + cable.R_i(2:end))/2;
cable.R_i_right(end) = inf; 

% Biophysics of cable
cable.V_rest = V_rest;

% Ion channels
cable.E_Na = E_Na;
cable.E_K  = E_K;
cable.E_L  = E_L;

cable.g_Na = g_Na;
cable.g_K  = g_K;
cable.g_L  = g_L;


switch modprmtr.mod
    case 1
        cable.z_ind_no_act = 1;                                     % Force activation term to zero at antidronic terminal
        z_ind_AP = 1;                                               % Location for AP to reach antidronic terminal for detection
        z_prop = abs(cable.z(N_comp) - cable.z(z_ind_AP));          % Propagation distance from AP initiation point (terminal) to AP detection point
    case 3
        cable.z_ind_no_act = 1;                             	% Disabling the terminals
        z_ind_AP = 1;
        % Distance for AP to pass; minimum of 10 mm or 2.5 times distance where
        % hyperpolarization is strongest (at sqrt(3/2) * H = 1.2 * H): 3*H
        z_prop = min(abs(median(cable.z) - cable.z(z_ind_AP))); % Assume AP initiation somewhere close to center of axon
end


%% E-field at cable

switch modprmtr.mod
    case 1
        if (modprmtr.ismicrofield == 1)
            E = 1;                      %  1 mV/cm
            Ex = 0 * ones_A;      % Tranverse E-field in mV/cm
            Ez = E * ones_A;      % Axial field in mV/cm
            d_phi_e = - Ez(1:end-1).*cable.dz(1:end-1)/2 - Ez(2:end).*cable.dz(2:end)/2;
        else
            ind_N = ceil(N/10);
            % ind_pos = mod(N, 10);
            E = E_z{ind_N};
            ind_pos = round( (1 - mod(N, 10)*0.02) * length(E));
            E = [E(1)*ones(1, length(E)-ind_pos), E(1:ind_pos)];
            Ex = 0 * ones_A;      % Tranverse E-field in mV/cm
            Ez = E;      % Axial field in mV/cm
            phi_e = cumsum(Ez*dz_vec);
            phi_e = -interp1(z_vec, phi_e, z);
            d_phi_e = diff(phi_e);    
        end
        

    case 3

    if (modprmtr.ismicrofield == 1)
        V = V1/f;
    else
        V = V2;
    end
    Ex = 0 * ones_A;      % Tranverse E-field in mV/cm
    phi_e = interp1(x1, V, z, 'linear', 0)/8;        
    d_phi_e = diff(phi_e);
end

stimulation.d_phi_e_left =  [ 0 ; d_phi_e ];    % Finite difference in potenial, in mV
stimulation.d_phi_e_right = [ d_phi_e ; 0 ];    % Finite difference in potenial, in mV

stimulation.ER_TP = kron(   (1+1./cable.TP_dim) .* abs( Ex ) .*...
                             cable.R , cos(theta)  );



%% Time and stimulation waveform

n_bf_start = 5;
t_start = - n_bf_start * dt;                                % Pulse on-set delay, in ms
t_end =   ceil( ( PW  + 5 + z_prop/v)/ dt ) * dt;           % Simulation end time, PW + ~ 5 ms for AP initiation + propagation time
t_vec = ( t_start : dt : t_end);                            % Time vector
if abs(round(PW/dt) - (PW/dt)) > 1e-3                       % If pulse width is not a multiple of dt
    t_vec = sort([t_vec,PW]);                                % Include PW in time vector
end    

stimulation.pulse_shape = zeros(size(t_vec));
stimulation.pulse_shape( (t_vec > 1e-6) & t_vec <= (PW + 1e-6) ) = 1;
stimulation.PW = PW;

%% Solver related parameters
solver.n_theta = N_theta;
solver.t_vec = t_vec;
solver.Temp = T;
solver.V_init = V_rest;
solver.h_func = str2func('simulate_cable_HH');

switch modprmtr.mod
    case 1
        E_rh = 100;          % mV/cm
        t_ch = 0.3;         % ms
        amp_init = E_rh / ( 1 - 2^ ( -PW/t_ch ) ) / (R_a/3e-4);        % mV/cm
        
        amp_init = amp_init * (1 + randn(1) * 0.02 );

        solver.plot_t_intv = 0.25:0.25:solver.t_vec(end);      % 250 us interval for plotting; 4 points per ms
        
        solver.thresh_find.unit_str = 'V/m';
        solver.thresh_find.unit_amp = 0.1; % 1 mV/cm = 0.1 V/m
        
        factor = sqrt(2);
  
    case 3
        % set initial search amplitude (based on RMG paper)
        t_ch = 0.36; % chronaxie from RMG; 102+-8 us
        k_rh = 1e5; % uA/cm^2; k is about 100 uA/mm^2 for 0.1 ms in RMG
        k_PW = k_rh / ( 1- 2^ (-PW/t_ch)); %uA/cm^2
        I_rh = 10;  %uA; I_0 is 25 uA for 0.1 ms in RMG
        I_0_PW = I_rh / ( 1- 2^ (-PW/t_ch)); %uA
        H = 30e-4;
        amp_init = - ( I_0_PW + k_PW * H^2 ) / (R_a/3e-4); % uA, cathodic current
        amp_init = amp_init * 0.2 * (1 + randn(1) * 0.02 );  % start at reduce amplitude and add some variation
        
        solver.plot_t_intv = 0.10:0.10:solver.t_vec(end);               % 100 us interval for plotting; 50 points per ms
                
        solver.thresh_find.unit_str = 'A';
        solver.thresh_find.unit_amp = 1e-6; % 1 uA = 1e-6 A
        factor = sqrt(2);
end

solver.thresh_find.amp_init = amp_init;
solver.thresh_find.amp_th_acc = 0.05e-2;             % 0.05%, accuracy of threshold finding
solver.thresh_find.factor = factor;
solver.thresh_find.range = 10^5; 
solver.thresh_find.phi_AP = 0;                      % mV; threshold definition, phi_m to cross 0 mV
solver.thresh_find.z_ind_AP = z_ind_AP;

end