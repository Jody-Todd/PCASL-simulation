% simulate effects of velocity
% calculate labeling efficiency

% Set Parameters

b0 = 3;

% pcasl parameters
rf_gap=1500e-06;                 % time between labeling pulses
label_dur=1.5;                    % labeling duration [s]
rf_fa=59;                       % flip angle of each label pulse [degrees]
rf_dur=500e-06;                 % duration of each labeling pulse
g_max=9;                        % max Gz during each labeling pulse [mT/m]
g_ave=1;                        % average gradient between rf pulses [mT/m]
labelPlaneThick=3e-03;          % thickness of labeling plane [m]

T1 = 1.55;                      % T1 blood [s]
T2 = 0.25;                      % T2 blood [s]
gamma = 42.58e06;               % gyromagnetic ratio [Hz/T]

% simulation parameters
v=20e-02;                       % velocity of blood spin [m/s]
dt=0.01e-03;                    % integration step size [s]  

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
Ngz = int32((gz.flatTime + gz.riseTime*2)/dt); % number of points in positive lobe
gz_wav = construct_g(gz.riseTime, gz.amplitude, Ngz, dt);
%figure; plot(gz_wav)

%%
%%%% tune gn_amp such that average gradient between rf pulses is g_ave %%%%
gn_amp = 127200; % [Hz/m]
gn_amp = 182000;
Ngn = floor((rf_gap - mr.calcDuration(gz))/dt); % number of points in negative lobe

% calculate average gradient between labeling pulses
[avg_grad, gn_wav] = get_avg_grad(gz_wav, gn_amp, Ngn, sys, dt);

fprintf('average gradient: %8f', avg_grad);
fprintf('\n g_ave: %8f \n', mr.convert(g_ave, 'mT/m'));
gz_full = [gz_wav gn_wav];

%%
figure; plot(gz_full);


% resample rf waveform to dt
rf_dt = rf.t(1):dt:rf.t(end);
rf_res = interp1(rf.t, rf.signal, rf_dt);

rf_padded = zeros(1,length(gz_full));
[gmax,idx] = max(gz_wav); % get index of max which corresponds to end of ramp period
rf_padded(idx-1:idx+length(rf_res)-2) = rf_res;

time_vec = 0:dt:dt*(length(gz_full)-1);

figure;
sgtitle('Single RF and Gz block', 'FontSize', 11)
subplot(211); plot(time_vec, rf_padded); title('RF'); ylabel('RF [Hz]');
subplot(212); plot(time_vec, gz_full); title('Gz'); ylabel('Gz [Hz/m]');
xlabel('Time [s]');


%seq.addBlock(rf, gz);
%seq.plot('timeDisp','s','showBlocks',1); %detailed view


% build rf and gradient waveforms such that the total labeling time equals
% the labeling duration
n = 1;
block_dur = length(gz_full)*dt;
total_dur = block_dur;
while total_dur <= label_dur
    total_dur = n*block_dur;
    n = n+1;
end

gz_block = repmat(gz_full, [1,n]);
rf_block = repmat(rf_padded, [1,n]);
block_time = time_vec(1):dt:(length(rf_block)-1)*dt + time_vec(1);

figure;
sgtitle('Full Labeling block', 'FontSize', 11)
subplot(211); plot(block_time, rf_block); title('rf'); ylabel('cycles/t');
subplot(212); plot(block_time, gz_block); title('gz'); ylabel('Hz/m');
xlabel('time [s]');



% loop through time points and peform a finite rotation

sim_time_vec = -2:dt:2;
n = length(sim_time_vec);
M_store = zeros(3,n);
M = [0; 0; 1];          
i = 1; % counter for labeling blocks
j = 1; % counter for M_store

for t = sim_time_vec
    if i == length(gz_full)
        i = 1;
    end
    g = gz_full(i);
    b1 = rf_padded(i);

    % Calculate Field Offsets
    freq_offset = v*g*t;  % Delta Omega [Hz]
    
    % In rotating frame: B_eff = [B1,  0,  freq_offset]. Assume B1
    % applied along x.
    Beff_vec = [mr.convert(b1,'Hz','T'); 0; freq_offset/gamma]; % [T]
    Omega_eff = norm(gamma*Beff_vec); % Magnitude of effective field
    
    % Perform a finite rotation using Rodrigues' Formula
    if Omega_eff > 0
        k = Beff_vec / norm(Beff_vec); % Normalized axis of rotation
        alpha = 2*pi*dt*Omega_eff; % Rotation angle [rad]
        
        % Rodrigues' Rotation Formula: v_rot = v cos(a) + (k x v) sin(a) + k(k.v)(1-cos(a))
        M = M*cos(alpha) + ...
            cross(k, M)*sin(alpha) + ...
            k*dot(k, M)*(1 - cos(alpha));
    end
    
    % Relaxation
    E2 = exp(-dt/T2);
    E1 = exp(-dt/T1);

    % perform T2 relaxation in x and y
    M(1) = M(1)*E2;
    M(2) = M(2)*E2;
    
    % peform T1 relaxation in z
    M0 = 1.0;
    M(3) = (M(3) - M0)*E1 + M0;
    
    % Store and Increment
    M_store(:,j) = M;
    j = j + 1;
    i = i + 1;
end

%%
% Visualization
figure;
plot(sim_time_vec, M_store(3, :), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Mz (Longitudinal Magnetization)');
title('Adiabatic Inversion Simulation');
grid on;
ylim([-1.1 1.1]);

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


function [avg_grad, gn_wav] = get_avg_grad(g_wav, gn_amp, Ngn, sys, dt)
% this function calculates the average gradient given one waveform and
% parameters to calculate a second waveform. this funciton can be used to
% determine the negative gradient lobe during rf_gap such that the average
% gradient between labeling pulses is g_ave.
% 
% Parameters:
%   g_wav: positive gradient waveform to include in the average.
%   gn_amp: amplitude of negative waveform
%   Ngn: number of negative waveform samples
%   sys: system parameters specified by mr.opts()
%   dt: sample time (can be different from gradRasterTime)

riseTime = abs(gn_amp)/sys.maxSlew;
riseTime = ceil(riseTime/dt)*dt;

% construct negative gradient waveform
gn_wav = construct_g(riseTime, -gn_amp, Ngn, dt);

% concatonate positive and negative waveforms and average
g_concat = [g_wav gn_wav];
avg_grad = mean(g_concat);

end

