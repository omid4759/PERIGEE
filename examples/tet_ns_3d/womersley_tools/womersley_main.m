% Rigid Womersley solution. Performs FFT on an inflow profile and computes
% the axial velocity as a function of r, t. 
% pressure:  p = A_n * exp(inwt)
% flow rate: Q = B_n * exp(inwt)


% Fluid properties
mu = 0.04;                           % dynamic viscosity (dyn * s / cm^2)
rho = 1.00;                          % density (g / cm^3)
nu = mu / rho;                       % kinematic viscosity (stokes)


% Flow rate (mL/s) measured in a pig's main pulmonary artery
T = 1.1; % period (s)
R = 0.3; % pipe radius (cm)
flow = [ 33.42, 56.19, 73.697, 96.721, 139.85, 164.46, 177.44, 196.25, ...
         198.77, 184.72, 162.09, 131.85, 91.057, 75.404, 62.991, 32.539, ...
         21.865, 28.182, 23.896, 19.457, 19.911, 13.432, 5.284, -1.0584 ];

% Number of Fourier modes (including the steady 0th mode)
n_modes = 2;

% Plot u_z(r, t)
B_n = womersley_velocity(flow / 50, nu, T, R, n_modes);

% Compute corresponding Fourier coefficients for pressure
A_n = fourier_modes_Q2p(B_n, rho, mu, omega, R, n_modes);

% Compute wall shear stress
wss = womersley_wss(A_n, nu, T, R, n_modes);
tawss = mean(wss);