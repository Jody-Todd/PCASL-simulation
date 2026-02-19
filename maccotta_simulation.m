% Adiabatic Inversion Simulation based on Maccotta et al., 1997

% Parameters
v = 20E-02;             % velocity of blood [m/s]
dt = 0.1E-03;           % time step [s]
b1 = 2.4E-06;           % RF field [T] (~24 mG)
g = 0.0025;             % gradient [T/m] (0.25 G/cm)
gamma = 42.58E06;       % gamma [Hz/T]
T1 = 1.0;               % T1 blood [s]
T2 = 0.2;               % T2 blood [s]

time_vec = -2:dt:2;     
n = length(time_vec);
M_store = zeros(3,n);

M = [0; 0; 1];          

i = 1;
for t = time_vec
    % Calculate Field Offsets
    freq_offset = gamma*v*g*t;  % Delta Omega [Hz]
    
    % In rotating frame: B_eff = [B1,  0,  freq_offset/gamma]. Assume B1
    % applied along x.
    Beff_vec = [b1; 0; freq_offset/gamma]; 
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
    M_store(:,i) = M;
    i = i + 1;
end

% Visualization
figure;
plot(time_vec, M_store(3, :), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Mz (Longitudinal Magnetization)');
title('Adiabatic Inversion Simulation');
grid on;
ylim([-1.1 1.1]);


