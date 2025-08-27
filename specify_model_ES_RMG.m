function [solver, stimulation, cable] = specify_model_ES_RMG(modprmtr)

PW_vec = [ kron([1e-2, 1e-1, 1e0], [1, 3.1]), 1e1];
% PW in ms, 10 us to 10 ms; >10 ms rheobase (steady state for linear
% membrane)

if strncmp(modprmtr.model_name,'UF',2)                  % Uniform Field
    mod_str = 'Axon terminal';
    load('E-field_profiles.mat', 'N', 'lambda', 'phi', 'A', 'dz', 'z_vec', 'E_z');
    dz_vec = dz;
    modprmtr.has_synapse = 0;
    L_straight = max(z_vec);         % Axon length, 10 mm in cm
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
    L_straight = 1;         % Axon length, 10 mm in cm

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
% Richardson, McIntyre & Grill 2000 (RMG)
T = 37;                     % Temperature, Degree Celcius

% Conductivities of intra- and extra-cellular spaces
sigma_i = 1/0.07;           % Intracellular conductivity, in mS/cm; resistivity: 70 Ohm*cm = 0.07 kOhm*cm
sigma_e = 1/0.5;            % Extracellular conductivity, in mS/cm; resistivity: 500 Ohm*cm = 0.5 kOhm*cm

if strncmp(modprmtr.model_size, '3um', 3)
    R_n = 1.65e-4;              % Nodal radius, in cm; 3.3 um diameter
    R_in = 3e-4;                % Internodal axon radius, in cm; 6 um fiber diameter (1 um diameter including myelin)
    N_lamella = 120;            % 120 myelin lamella; 2 membranes per lamella
    L_in = 1150e-4;             % Internodal length, in cm; 1150 um
    N_in = 1150;             	% 1150 compartments per internode,
elseif strncmp(modprmtr.model_size, '0.3um', 5)
    R_n = 0.165e-4;             % Nodal radius, in cm; 0.33 um diameter
    R_in = 0.3e-4;              % Internodal axon radius, in cm; 0.6 um fiber diameter (1 um diameter including myelin)
    N_lamella = 12;             % 12 myelin lamella; 2 membranes per lamella
    L_in = 115e-4;              % Internodal length, in cm; 115 um
    N_in = 115;             	% 115 compartments per internode,
end


% Nodal membrane parameters
c_n = 2;                    % Nodal membrane capacitance, in uF/cm^2
l_n = 1e-4;                 % Nodal length, in cm; 1 um
% "Cellular" time constant for cylindrical nodal segment
% tau_c = R_n * c_n * (sigma_i^-1 + sigma_e^-1) = 188.1 ns
% 3 * tau_c = 0.56 us  ->  reaching 95.0% of steady state TP
% 5 * tau_c = 0.94 us  ->  reaching 99.3% of steady state TP

% Ion channel parameters for nodes, in mV & mS/cm^2;
V_rest = -82;                       % in MRG 2002: -80;
E_Na = 50;        g_Na = 3000;      g_Nap = 5;  % fast & persistent Na
E_K  = -84;       g_K  = 80;        % Slow K
E_L  = -83.38;    g_L  = 80;        % Leakage
% g_bar = 104.3 mS/cm^2 at rest -> r_m = 0.0096 kOhm*cm^2
% tau_n = 19.2 us

% Internodal membrane parameters
c_in = 1/2/N_lamella;       % Internodal total myelin capacitance, in uF/cm^2
g_in = 1/2/N_lamella;       % Internodal total myelin conductance, in mS*cm^2
l_in = L_in / N_in;         % Internodal compartment length 1, um, in cm

% lambda:
[m,~,h,~,p,~,s,~] = RMG(V_rest, T);
g_bar = ( g_Na .* m.^3 .* h + g_Nap .* p.^3 )  + (g_K.* s) + g_L;

lambda_m = sqrt( (sigma_i * R_in) / (2*g_in ) );
tau_m = c_in / g_in;
lambda_n = sqrt( (sigma_i * R_n ) / (2*g_bar) );
tau_n = c_n / g_bar;

lambda_DC =  1 / sqrt( (1 - l_n / (l_n+L_in)) / lambda_m^2 + ...
    l_n / (l_n+L_in)  / lambda_n^2 );

tau = lambda_DC^2 *  ( (1 - l_n / (l_n+L_in)) / lambda_m^2 * tau_m + ...
    l_n / (l_n+L_in)  / lambda_n^2 * tau_n );

% Emperical parameter from simulation/RMG paper
v = 0.61;                    % cm/ms; conduction speed 61 m/s = 6.1 cm/ms from RMG paper
v = v * (R_in/3e-4);

%% Specify cable
N_theta = 15;               % Discretization of azimuthal angle
d_theta =  pi / N_theta;    % Interval for integration, in radian
theta = linspace( d_theta/2, pi - d_theta/2, N_theta );         % Integration points between 0 and pi, row vector

switch modprmtr.mod
    case 1
        N_node = floor(L_straight/(L_in+l_n)) + 1;              % Number of nodal compartments
        N_comp = N_node + (N_node - 1) * N_in;       % Number of all compartments
    case 3

        N_node = floor(L_straight/(L_in+l_n)) + 1;              % Number of nodal compartments
        N_comp = N_node + (N_node - 1) * N_in;       % Number of all compartments

        L_straight = (N_node-1) * (L_in+l_n) + l_n;
        axon_pos_length = max(x1)*(N-1)/40;
        axon_neg_length = L_straight - axon_pos_length;  % +2 mm

        z = (-axon_neg_length :dz: axon_pos_length )';     % axial coordinates, in cm. Row vector

end

switch modprmtr.mod
    case 1

        N_myl_comp = N_in * 10;
        N_comp = N_comp + N_myl_comp * 1;

    case 3
        N_myl_comp = N_in * 10;
        N_comp = N_comp + N_myl_comp * 1;
end

Empty_N = zeros(N_node,1);              % Column vector; empty array
Empty_th = zeros(N_node,N_theta);       % Column vector expanded; empty array
node = struct(  'cable_ID',Empty_N,...  % the index of the node within the entire cable
    'TP_dim',Empty_N,...    % Dimension: 1 for cylindrical; 2 for sphercial
    'TP_weight',Empty_th... % Integration weights, row vector for each node
    );
Empty_C = zeros(N_comp, 1);             % Column vector; empty array
cable = struct( 'N_comp',N_comp,...     % Number of compartments
    'z',Empty_C,...         % Center coordinates of compartment, cm (along local longitudinal axis)
    'dz',Empty_C,...        % Compartment length, cm
    'R',Empty_C,...         % Compartment radius, cm
    'R_i',Empty_C,...       % Compartment axial resistance, kOhm
    'R_i_left',Empty_C,...  % Axial resistance to left neighbor, kOhm
    'R_i_right',Empty_C,... % Axial resistance to right neighbor, kOhm
    'Area',Empty_C,...      % Compartment Area, cm^2
    'C_m',Empty_C,...       % Compartment capacitance, uF
    'is_node',Empty_C,...   % Is compartment a node? logical 0/1
    'N_nodes',N_node,...    % Number of nodes
    'node',node...          % Nodal data
    );


z_in_rel = l_in * (-0.5 + (1 : N_myl_comp));
ind_in = (1:N_myl_comp) ;
cable.z(ind_in) = z_in_rel;
cable.dz(ind_in) = l_in;                                % Compartments' length, in cm
cable.R(ind_in) = R_in;                                 % Compartments' radius, in cm
cable.R_i(ind_in) = l_in / (sigma_i * pi * R_in^2);     % Compartments' axial resistance, in kOhm
cable.Area(ind_in) = 2 * pi * R_in * l_in;              % Compartments' area; cm^2
cable.C_m(ind_in) = c_in * 2 * pi * R_in * l_in;        % Compartments' membrane capacitance, in uF
cable.is_node(ind_in) = false;                          % Compartments are internode

z_in_rel = l_in * ( (1:N_in) - 0.5 );
% Relavtive axial coordinates of the center of N_in internodal compartments
% within one internode from its left end

for ii = 1 : N_node      % Add node and internode
    % 1 nodal compartments
    ind_n = N_myl_comp + (1 + N_in) * (ii-1) + 1;           % Index of compartment in the cable: (ii-1) nodes and internodes to the left
    cable.z(ind_n) = N_myl_comp * l_in + (l_n + L_in) * (ii-1) + l_n/2;         % Center coordinates of compartment, cm
    cable.dz(ind_n) = l_n;                                  % Compartment's length, in cm
    cable.R(ind_n) = R_n;                                   % Compartment's radius, in cm
    cable.R_i(ind_n) = l_n / (sigma_i * pi * R_n^2);        % Compartment's axial resistance, in kOhm
    cable.Area(ind_n) = 2 * pi * R_n * l_n;                 % Compartment's area; cm^2
    cable.C_m(ind_n) = c_n * 2 * pi * R_n * l_n;            % Compartment's membrane capacitance, in uF
    cable.is_node(ind_n) = true;                            % Compartment is node
    cable.node.cable_ID(ii) = ind_n;                    % The index of the node within the entire cable
    cable.node.TP_dim(ii) = 1;                          % Dimension parameter for cylindrical compartments is 1
    cable.node.TP_weight(ii,:) = d_theta / pi;          % Each patch of the membrane is only a fraction of the total area
    if ii == N_node
        break                                               % No internode after last nodal compartment
    end

    % N_in internodal compartments
    ind_in = N_myl_comp + ii + N_in * (ii-1) + (1:N_in);	% Index of compartments in the cable: ii nodes and  (ii-1) internodes to the left
    cable.z(ind_in) = N_myl_comp * l_in + l_n * ii + L_in * (ii-1) + z_in_rel;  % Center coordinates of compartments, cm
    cable.dz(ind_in) = l_in;                                % Compartments' length, in cm
    cable.R(ind_in) = R_in;                                 % Compartments' radius, in cm
    cable.R_i(ind_in) = l_in / (sigma_i * pi * R_in^2);     % Compartments' axial resistance, in kOhm
    cable.Area(ind_in) = 2 * pi * R_in * l_in;              % Compartments' area; cm^2
    cable.C_m(ind_in) = c_in * 2 * pi * R_in * l_in;        % Compartments' membrane capacitance, in uF
    cable.is_node(ind_in) = false;                          % Compartments are internode
end

switch modprmtr.mod
    case 1

        cable.z_ind_no_act = 1;                                     % Force activation term to zero at antidronic terminal
        z_ind_AP = 1 + N_myl_comp + (1 + N_in) * 0;                 % Location for AP to reach antidronic terminal for detection
        z_prop = abs(cable.z(N_comp) - cable.z(z_ind_AP));          % Propagation distance from AP initiation point (terminal) to AP detection point
        cable.z = cable.z - max(cable.z-max(z_vec));

    case 3
        cable.z_ind_no_act = 1;                                     % Force activation term to zero at antidronic terminal
        z_ind_AP = 1 + N_myl_comp + (1 + N_in) * 0;                 % Location for AP to reach antidronic terminal for detection
        z_prop = abs(cable.z(N_comp) - cable.z(z_ind_AP));          % Propagation distance from AP initiation point (terminal) to AP detection point
        cable.z = cable.z - (max(cable.z)-max(z));
end


% Axial resistance; sealed ends are reflected by d_phi_e_left/right = 0 at terminals
cable.R_i_left(1) = inf;
cable.R_i_left(2:end)    = ( cable.R_i(1:end-1) + cable.R_i(2:end))/2;
cable.R_i_right(1:end-1) = ( cable.R_i(1:end-1) + cable.R_i(2:end))/2;
cable.R_i_right(end) = inf;

% Biophysics of cable
cable.V_rest = V_rest;
cable.g_in = g_in;

% Nodal ion channels
cable.node.E_Na = E_Na;
cable.node.E_K  = E_K;
cable.node.E_L  = E_L;

cable.node.g_Na = g_Na;
cable.node.g_Nap = g_Nap;
cable.node.g_K  = g_K;
cable.node.g_L  = g_L;

%% E-field at cable

switch modprmtr.mod
    case 1
        if (modprmtr.ismicrofield == 1)
            E = 1;                      %  1 mV/cm
            Ex = 0 * Empty_C;      % Tranverse E-field in mV/cm
            Ez = E * ones(size(Empty_C));      % Axial field in mV/cm
            d_phi_e = - Ez(1:end-1).*cable.dz(1:end-1)/2 - Ez(2:end).*cable.dz(2:end)/2;
        else
            ind_N = ceil(N/10);
            % ind_pos = mod(N, 10);
            E = E_z{ind_N};
            ind_pos = round( (1 - mod(N, 10)*0.02) * length(E));
            E = [E(1)*ones(1, length(E)-ind_pos), E(1:ind_pos)];
            Ex = 0 * Empty_C;      % Tranverse E-field in mV/cm
            Ez = E;      % Axial field in mV/cm
            phi_e = cumsum(Ez*dz_vec);
            phi_e = -interp1(z_vec, phi_e, cable.z, 'linear', 'extrap');
            d_phi_e = diff(phi_e);
            d_phi_e(cable.z<0) = -1 * median(diff(cable.z));
        end

    case 3
        I = 1;                                      % 1 uA; scaling performed by waveform

        if (modprmtr.ismicrofield == 1)
            V = V1/f;
        else
            V = V2;
        end
        Ex = 0 * Empty_C;      % Tranverse E-field in mV/cm
        phi_e = interp1(x1, V, cable.z, 'linear', V(1))/8;
        d_phi_e = diff(phi_e);
        Ez = Ex;
end

stimulation.d_phi_e_left =  [ 0 ; d_phi_e ];    % Finite difference in potenial, in mV
stimulation.d_phi_e_right = [ d_phi_e ; 0 ];    % Finite difference in potenial, in mV

stimulation.ER_TP = kron(   (1+1./cable.node.TP_dim) .* abs( Ex(cable.node.cable_ID) ) .*...
    cable.R(cable.node.cable_ID) , cos(theta)  );

%% Time and stimulation waveform
dt = 1e-3;              % Time step 1 us, in ms

n_bf_start = 5;
t_start = - n_bf_start * dt;                                % Pulse on-set delay, in ms
t_end =   ceil( ( PW  + 1 + z_prop/v)/ dt ) * dt;           % Simulation end time, PW + ~ 1 ms for AP initiation + propagation time
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
solver.h_func = str2func('simulate_cable_RMG');

switch modprmtr.mod
    case {1,2}
        E_rh = 5;           % mV/cm
        t_ch = 0.3;         % ms
        amp_init = E_rh / ( 1 - 2^ (-PW/t_ch)) / (R_n / 1.65e-4);        % mV/cm

        amp_init =  amp_init * (1 + randn(1) * 0.02 );

        solver.plot_t_intv = 0.1:0.1:solver.t_vec(end);      % 50 us interval for plotting; 20 points per ms

        solver.thresh_find.unit_str = 'V/m';
        solver.thresh_find.unit_amp = 0.1; % 1 mV/cm = 0.1 V/m
        factor = (2)^(1/2);
    case 3
        % set initial search amplitude (based on RMG paper)
        t_ch = 0.1;     % chronaxie from RMG; 102+-8 us
        k_rh = 2000;	% uA/cm^2; k is about 100 uA/mm^2 for 0.1 ms in RMG
        k_PW = k_rh / ( 1 - 2^ (-PW/t_ch));     %uA/cm^2
        I_rh = 2;       %uA; I_0 is 25 uA for 0.1 ms in RMG
        I_0_PW = I_rh / ( 1 - 2^ (-PW/t_ch));   %uA

        amp_init = - ( I_0_PW + k_PW * 30e-4^2 ) / (R_n / 1.65e-4);  % uA, cathodic current
        amp_init = amp_init * 0.2 *(1 + randn(1) * 0.02 );  % start at reduce amplitude and add some variation

        % if H <= 0.5
        solver.plot_t_intv = 0.05:0.05:solver.t_vec(end);             % 20 us interval for plotting; 50 points per ms
        % else
        %     solver.plot_t_intv = [0.05:0.05:0.9,1:0.1:solver.t_vec(end)];	% 100 us interval for plotting; 10 points per ms
        % end

        solver.thresh_find.unit_str = 'A';
        solver.thresh_find.unit_amp = 1e-6; % 1 uA = 1e-6 A
        factor = sqrt(2);
end

solver.thresh_find.amp_init = amp_init;
solver.thresh_find.amp_th_acc = 0.05e-2;             % 0.5%, accuracy of threshold finding
solver.thresh_find.factor = factor;
solver.thresh_find.range = 10^4;
solver.thresh_find.phi_AP = 0;                      % mV; threshold definition, phi_m to cross 0 mV
solver.thresh_find.z_ind_AP = z_ind_AP;

end
