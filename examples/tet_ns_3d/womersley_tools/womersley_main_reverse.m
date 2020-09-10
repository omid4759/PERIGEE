% Rigid Womersley solution. Start with Fourier coefficients for pressure
% and reconstruct the flow rate.

% pressure:  p = A_n * exp(inwt)
% flow rate: Q = B_n * exp(inwt)

% Fluid properties
mu = 0.04;                           % dynamic viscosity (dyn * s / cm^2)
rho = 1.00;                          % density (g / cm^3)
nu = mu / rho;                       % kinematic viscosity (stokes)

T = 1.1; % period (s)
R = 0.3; % pipe radius (cm)

% base angular frequency
omega = 2 * pi / T;

% Time vector
N = 24;                              % number of samples per period
t = linspace(0, T, N);

% Fourier coefficients for pressure (including the steady 0th mode)
A_n = [-21.046886507506900, -33.010163920166285 + 42.933217768314500 * 1i];
% A_n = [10, 10 + 10 * 1i];

% Number of Fourier modes (including the steady 0th mode)
n_modes = length(A_n);

% Compute flow rate Fourier coefficients
B_n = fourier_modes_p2Q(A_n, rho, mu, omega, R, t, n_modes);

% Compute and plot v_z(r, t)
womersley_velocity_reverse(B_n, nu, T, R, omega, t, n_modes);

% Compute wall shear stress
wss = womersley_wss(A_n, nu, T, R, n_modes);
tawss = mean(wss);
