
function womersley_velocity_reverse(B_n, nu, T, R, omega, t, n_modes)

% ---------- Reconstruct the flow rate ----------
flow_recon = zeros(1, length(t));
for n = 1 : n_modes
    coef = (n - 1) * omega;
    flow_recon = flow_recon + B_n(n) * (cos(coef * t) + 1i * sin(coef * t));
end
figure; plot(t, real(flow_recon));                          % reconstructed flow
xlim([0, T])
xlabel('Time (s)'); ylabel('Reconstructed Flow (mL/s)')
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

