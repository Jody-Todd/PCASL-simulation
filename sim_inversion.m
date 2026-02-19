function [sim_time_vec, M_store] = sim_inversion(gz_shape, rf_shape, params)
% Simulates the adiabatic inversion of a spin flowing through the labeling 
% plane during PCASL labeling. This funciton performs a numerical 
% integration of the bloch equations by performing finite rotations around
% the effective magnetic field Beff in the rotating reference frame 
% (roational frequency = omega_rf).
% 
% Parameters:
%   gz_shape: z-gradient waveform applied during a single RF pulse in Hz/m.
%   rf_shape: single rf pulse in T. expected to be same size and sampling time 
%       as gz_shape.
%   params: structure with fields T1 (longitudinal relaxation time of blood
%       [s]), T2 (transverse relaxation time of blood [s]), gamma (gyromagnetic
%       ratio [Hz/T]), v (velocity [m/s]), dt (step time [s]), sim_time (postive
%       number that indicates the end of the simulation time. simulation 
%       will run from -sim_time:dt:sim_time).
%
% Returns:
%   sim_time_vec: simulation time vector.
%   M_store: evolution of magnetization [Mx; My; Mz] as it experiences adiabatic
%   inversion.
%

dt=params.dt;
E2 = exp(-dt/params.T2);
E1 = exp(-dt/params.T1);
gamma=params.gamma;
v=params.v;
t1=-params.sim_time;
t2=params.sim_time;

sim_time_vec = t1:dt:t2;
M_store = zeros(3,length(sim_time_vec));
M = [0; 0; 1];
M0 = 1.0;
i = 1; % counter for labeling blocks
j = 1; % counter for M_store

for t = sim_time_vec
    if i == length(gz_shape)
        i = 1;
    end

    % Calculate Field Offsets
    freq_offset = v*gz_shape(i)*t;  % Delta Omega [Hz]
    
    % In rotating frame: B_eff = [B1,  0,  freq_offset]. Assume B1
    % applied along x.
    Beff_vec = [rf_shape(i); 0; freq_offset/gamma]; % [T]
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

    % perform T2 relaxation in x and y
    M(1) = M(1)*E2;
    M(2) = M(2)*E2;
    
    % peform T1 relaxation in z
    M(3) = (M(3) - M0)*E1 + M0;
    
    % Store and Increment
    M_store(:,j) = M;
    j = j + 1;
    i = i + 1;
end



