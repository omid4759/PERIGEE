% Ingrid Lan
% ME 285 HW4

function B_n = womersley_velocity(flow, nu, T, R, n_modes)

% base angular frequency
omega = 2 * pi / T;
assignin('base', 'omega', omega);

% Time vector
N = floor(length(flow));
t = linspace(0, T, N);

% Discrete fast Fourier transform
Q_n = fft(flow);

% Compute the Fourier coefficients B_n for non-negative frequencies
B_n = Q_n / N;
B_n = B_n(1 : N / 2 + 1);
B_n(2 : end - 1) = 2 * B_n(2 : end - 1);

% Truncate to first n modes (including the steady 0th mode)
B_n = B_n(1 : n_modes);

% ------- Plot the input vs. reconstructed flow --------
figure; hold on; plot(t, flow);         % input flow
flow_recon = zeros(1, length(t));
for n = 1 : n_modes
    coef = (n - 1) * omega;
    flow_recon = flow_recon + B_n(n) * (cos(coef * t) + 1i * sin(coef * t));
end
plot(t, real(flow_recon));              % reconstructed flow
xlabel('Time (s)'); ylabel('Flow (mL/s)')
legend('Input', 'Reconstructed');
set(gca, 'FontSize', 12)

% ------- Plot axial velocity as a function of r, t -------
r = (-R : 0.01 : R)';                % radius vector
t_ds = 0 : T / 8 : T;                % downsampled time vector

v_z = zeros(length(r), length(t_ds));
figure; hold on;
legend_list = cell(1, length(t_ds));
line_style = {'k-', 'k--', 'k.-', 'r-', 'r--', 'r.-', 'b-', 'b--', 'b.-'};

for ii = 1 : length(t_ds)
    
    % Initialize with Poiseuille solution (0th mode)
    v_z(:, ii) = 2 * B_n(1) * (R^2 - r.^2) / (pi * R^4);
    
    for k = 2 : n_modes
        n = k - 1;
        alpha_n = R * sqrt(n * omega / nu);
        gamma_n = alpha_n * 1j^(3 / 2);
        v_z(:, ii) = v_z(:, ii) + B_n(k) * gamma_n  / (pi * R^2) * ...
              (besselj(0, gamma_n) - besselj(0, gamma_n * r / R)) / ...
              (besselj(0, gamma_n) * gamma_n - 2 * besselj(1, gamma_n)) * ...
              exp(1j * n * omega * t_ds(ii));
    end
    plot(real(v_z(:, ii)), r / R, line_style{ii});
    legend_list{ii} = ['t = ', num2str(t_ds(ii)), ' s'];
end
legend(legend_list, 'Location', 'Best');
xlabel('u_z (cm / s)'); ylabel('r / R');
set(gca, 'FontSize', 12)

% EOF