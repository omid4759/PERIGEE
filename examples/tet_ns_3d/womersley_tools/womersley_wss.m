
function wss = womersley_wss(A_n, nu, T, R, n_modes)

% base angular frequency
omega = 2 * pi / T;

t_ds = 0 : T / 10 : T;                % downsampled time vector

wss = -A_n(1) * R / 2;

for k = 2 : n_modes
    n = k - 1;
    alpha_n = R * sqrt(n * omega / nu);
    gamma_n = alpha_n * 1j^(3 / 2);
    
    wss = wss + A_n(k) * R * 1j^(1/2) / alpha_n * ...
          besselj(1, gamma_n) / besselj(0, gamma_n) * exp(1j*n*omega*t_ds);
      
    wss = real(wss);
end
