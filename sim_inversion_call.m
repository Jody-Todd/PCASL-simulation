%% Set Parameters

b0 = 3;

% pcasl parameters
rf_gap=1500e-06;                 % time between labeling pulses
label_dur=1.5;                    % labeling duration [s]
rf_fa=40;                       % flip angle of each label pulse [degrees]
rf_dur=500e-06;                 % duration of each labeling pulse
g_max=9;                        % max Gz during each labeling pulse [mT/m]
g_ave=1;                        % average gradient between rf pulses [mT/m]
labelPlaneThick=3e-03;          % thickness of labeling plane [m]

params.T1 = 1.55;                      % T1 blood [s]
params.T2 = 0.25;                      % T2 blood [s]
params.gamma = 42.58e06;               % gyromagnetic ratio [Hz/T]

% simulation parameters
params.v=20e-02;                       % velocity of blood spin [m/s]
params.dt=0.01e-03;                    % integration step size [s]  
params.sim_time=2;                     % simulation will occur over -sim_time:dt:sim_time [s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set system limits
sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
    'MaxSlew', 200, 'SlewUnit', 'T/m/s', ... 
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, 'B0', b0);

seq=mr.Sequence(sys);           % Create a new sequence object

% create a single labeling pulse and slice select gradient
N = round(rf_dur/sys.rfRasterTime);
bw = mr.convert(g_max, 'mT/m', 'Hz/m')*labelPlaneThick; % [Hz]
[rf, gz] = mr.makeArbitraryRf(hann(N), rf_fa*pi/180, 'sliceThickness', labelPlaneThick, 'bandwidth', bw, 'use', 'preparation', 'system', sys);

% construct gz waveform
Ngz = int32((gz.flatTime + gz.riseTime*2)/params.dt); % number of points in positive lobe
gz_wav = construct_g(gz.riseTime, gz.amplitude, Ngz, params.dt);

%% tune gn_amp such that average gradient between rf pulses is g_ave %%
gn_amp = 186680; % [Hz/m]
Ngn = floor((rf_gap - mr.calcDuration(gz))/params.dt); % number of points in negative lobe

% construct negative gradient lobe of given amplitude
gn_wav = construct_g(gz.riseTime, -gn_amp, Ngn, params.dt);

g_concat = [gz_wav gn_wav];
avg_grad = mean(g_concat);

if avg_grad > mr.convert(g_ave, 'mT/m')
    fprintf('\ngn amplitude too small\n')
elseif avg_grad < mr.convert(g_ave, 'mT/m')
    fprintf('\ngn amplitude too big\n')
end
fprintf('average gradient: %8f', avg_grad);
fprintf('\ng_ave: %8f \n', mr.convert(g_ave, 'mT/m'));
%%

gz_full = [gz_wav gn_wav];
figure; plot(gz_full);

% resample rf waveform to dt
rf_dt = rf.t(1):params.dt:rf.t(end);
rf_res = interp1(rf.t, rf.signal, rf_dt);

rf_padded = zeros(1,length(gz_full));
[gmax,idx] = max(gz_wav); % get index of max which corresponds to end of ramp period
rf_padded(idx-1:idx+length(rf_res)-2) = rf_res;
rf_padded = mr.convert(rf_padded, 'Hz', 'T');

time_vec = 0:params.dt:params.dt*(length(gz_full)-1);

figure;
sgtitle('Single RF and Gz block', 'FontSize', 11)
subplot(211); plot(time_vec, rf_padded); title('RF'); ylabel('RF [T]');
subplot(212); plot(time_vec, gz_full); title('Gz'); ylabel('Gz [Hz/m]');
xlabel('Time [s]');


% build rf and gradient waveforms such that the total labeling time equals
% the labeling duration
n = 1;
block_dur = length(gz_full)*params.dt;
total_dur = block_dur;
while total_dur <= label_dur
    total_dur = n*block_dur;
    n = n+1;
end

gz_block = repmat(gz_full, [1,n]);
rf_block = repmat(rf_padded, [1,n]);
block_time = time_vec(1):params.dt:(length(rf_block)-1)*params.dt + time_vec(1);

figure;
sgtitle('Full Labeling block', 'FontSize', 11)
subplot(211); plot(block_time, rf_block); title('RF'); ylabel('RF [T]');
subplot(212); plot(block_time, gz_block); title('Gz'); ylabel('Gz [Hz/m]');
xlabel('time [s]');


%% simulate adiabatic inversion for a few different velocities

[sim_time_vec, M_v20] = sim_inversion(gz_full, rf_padded, params);
params.v=10e-02;
[sim_time_vec, M_v10] = sim_inversion(gz_full, rf_padded, params);
params.v=5e-02;
[sim_time_vec, M_v5] = sim_inversion(gz_full, rf_padded, params);
params.v=30e-02;
[sim_time_vec, M_v30] = sim_inversion(gz_full, rf_padded, params);
params.v=50e-02;
[sim_time_vec, M_v50] = sim_inversion(gz_full, rf_padded, params);
params.v=100e-02;
[sim_time_vec, M_v100] = sim_inversion(gz_full, rf_padded, params);
params.v=0;
[sim_time_vec, M_v0] = sim_inversion(gz_full, rf_padded, params);

% Visualization
figure;
plot(sim_time_vec, M_v20(3, :), 'LineWidth', 1); hold on;
plot(sim_time_vec, M_v30(3, :), 'LineWidth', 1); hold on;
plot(sim_time_vec, M_v50(3, :), 'LineWidth', 1); hold on;
plot(sim_time_vec, M_v5(3, :), 'LineWidth', 1); hold on;
plot(sim_time_vec, M_v100(3, :), 'LineWidth', 1); hold on;
plot(sim_time_vec, M_v0(3, :), 'LineWidth', 1); hold on;

legend('20 cm/s', '30 cm/s', '50 cm/s', '5 cm/s', '100 cm/s', '0 cm/s')

xlabel('Time (s)');
ylabel('Mz (Longitudinal Magnetization)');
title('Adiabatic Inversion Simulation');
grid on;
ylim([-1.1 1.1]);


%% Efficiency Calculation for single velocity
eff = get_spin_efficiency(0.2, gz_full, rf_padded, params);

%% Efficiency calculation for a vessel with a velocity distribution

max_velocities = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1]; % maximum velocity at the center of the vessel [m/s]
efficiencies = zeros(1, length(max_velocities));
params.sim_time = 2;

% velocity-weighted efficieny integral
integrand = @(v) v .* get_spin_efficiency(v, gz_full, rf_padded, params);

i = 1;
for v_max = max_velocities

    % 'ArrayValued', true is required because get_spin_efficiency only takes scalar 'v' inputs
    integral_result = integral(integrand, 0.0001, v_max, 'ArrayValued', true, 'AbsTol', 0.00001);
    
    % Normalize by the laminar flow weighting factor 
    efficiencies(i) = 2*integral_result / (v_max^2);
    
    % Display the result
    fprintf('Vessel Labeling Efficiency for max velocity %.2f: %.2f%%\n', v_max, efficiencies(i) * 100);

    i = i + 1;
end

%% Visualization 
figure;
plot(max_velocities*100, efficiencies, 'LineWidth', 1);
xlabel('Maximum velocity in a vessel [cm/s]');
ylabel('Efficiency [%]');
ylim([0,1]);


%% Plot x and y magnetization

figure;
plot(sim_time_vec, M_store(1, :), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Mx (Transverse Magnetization)');
title('Adiabatic Inversion Simulation');
grid on;
ylim([-1.1 1.1]);

figure;
plot(sim_time_vec, M_store(2, :), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('My (Transverse Magnetization)');
title('Adiabatic Inversion Simulation');
grid on;
ylim([-1.1 1.1]);



function eff = get_spin_efficiency(v, gz_shape, rf_shape, params)
    % wrapper function for sim_inversion with specific velocity
    if v == 0
        eff = 0;
        return;
    end
    
    M0 = 1;

    % update the velocity parameter for this specific run
    params.v = v;
    
    % run Bloch simulation
    [~, M_store] = sim_inversion(gz_shape, rf_shape, params);
    
    % correct for T1 decay after the spin crosses the labeling plane (t=0)     
    T1_decay_factor = exp(-params.sim_time / params.T1);

    % extract the final longitudinal magnetization
    Mz_final = M_store(3, end);

    % calculate normalized efficiency (Perfect inversion = 1.0)
    % We negate Mz_final because a successful inversion results in negative magnetization.
    %eff = -(M0 + (Mz_final - M0)/T1_decay_factor); 
    %eff = (M0 - Mz_final)/(2*M0*T1_decay_factor);
    Mz_corrected = M0 + (Mz_final - M0)/T1_decay_factor;
    eff = (M0 - Mz_corrected)/(2*M0);
end



function g_wav = construct_g(riseTime, amp, N, dt)
% construct a trapezoidal waveform given riseTime, amplitude (amp), total
% number of samples (N), and sample time dt. instead of dt, one could use
% sys.gradRasterTime.

g_wav = zeros(1,N);

n_rise = floor(riseTime/dt); % number of rise samples
rise_samples = double(0:n_rise-1)*dt; % timing for rise samples
rise_slope = amp/(riseTime); % slope

n_flat = N - 2*n_rise; % number of flat samples

fall_samples = double(n_rise+n_flat:length(g_wav))*dt; % timing for fall samples
b = amp + rise_slope*fall_samples(1); % intercept

g_wav(1:n_rise) = rise_samples*rise_slope; % rise
g_wav(n_rise+1:n_rise+n_flat-1) = amp; % flat
g_wav(n_rise+n_flat:end) = fall_samples*-rise_slope + b; % fall

end
