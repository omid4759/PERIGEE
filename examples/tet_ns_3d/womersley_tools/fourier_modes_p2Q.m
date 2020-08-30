% Compute flow rate Fourier coefficients B_n from pressure Fourier 
% coefficients A_n

function B_n = fourier_modes_p2Q(A_n, rho, mu, omega, R, t, n_modes)
% p = A_n * exp(inwt)
% Q = B_n * exp(inwt)

B_n = zeros(1, n_modes);
B_n(1) = -1.0 * A_n(1) * pi* R^4 / (8.0 * mu);      % steady 0th mode

for k = 2 : n_modes
    n = k - 1;
    Omega = sqrt(rho * n * omega / mu) * R;         % Womersley number
    Gamma = 1i^1.5 * Omega;
    
    B_n(k) = -A_n(k) * pi * R^2 / (1i * rho * n * omega) * ...
             ( 1 - 2.0 * besselj(1, Gamma) / (Gamma * besselj(0, Gamma)) );
end
