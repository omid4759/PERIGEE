% Compute pressure Fourier coefficients A_n from flow rate Fourier 
% coefficients B_n

function A_n = fourier_modes_Q2p(B_n, rho, mu, omega, R, n_modes)
% p = A_n * exp(inwt)
% Q = B_n * exp(inwt)

A_n = zeros(1, n_modes);
A_n(1) = -1.0 * B_n(1) * 8.0 * mu / (pi* R^4);      % steady 0th mode

for k = 2 : n_modes
    n = k - 1;
    Omega = sqrt(rho * n * omega / mu) * R;         % Womersley number
    Gamma = 1i^1.5 * Omega;
    
    A_n(k) = B_n(k) * mu * Omega^2 / ...
             ( 1i * pi * R^4 * (1.0 - 2.0 * besselj(1, Gamma) / ...
               (Gamma * besselj(0, Gamma)) ));
end

% EOF