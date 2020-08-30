% manufactured solution for incompressible Navier-Stokes 3D
clear all; clc;

syms x y z t mu rho0 R k omega;

% Given analytic manufactured solution
u = 0.0;

v = 0.0;

Omega = sqrt(rho0 * omega/mu) * R;
coef = 1i * k / (rho0 * omega);
xi = (1i^1.5)*Omega*sqrt(x*x+y*y) / R;

w = coef * (1.0 - besselj(0,xi)/besselj(0,1i^1.5*Omega))*exp(1i*omega*t); %k*(x*x+y*y - R*R)/(4*mu);

p = k * z * exp(1i*omega*t);

% Compute first order derivatives
u_x = diff(u,x); u_y = diff(u,y); u_z = diff(u,z);
v_x = diff(v,x); v_y = diff(v,y); v_z = diff(v,z);
w_x = diff(w,x); w_y = diff(w,y); w_z = diff(w,z);

u_t = diff(u,t); v_t = diff(v,t); w_t = diff(w,t);

u_xx = diff(u_x, x); u_yy = diff(u_y, y); u_zz = diff(u_z, z);
v_xx = diff(v_x, x); v_yy = diff(v_y, y); v_zz = diff(v_z, z);
w_xx = diff(w_x, x); w_yy = diff(w_y, y); w_zz = diff(w_z, z);

p_t = diff(p, t); p_x = diff(p, x); p_y = diff(p, y); p_z = diff(p, z);

p_t = simplify(p_t);
p_x = simplify(p_x);
p_y = simplify(p_y);
p_z = simplify(p_z);

div_velo = u_x + v_y + w_z;

% Print divergence of velocity on screen
div_velo = simplify(div_velo)

% Compute grad dot velocity
u_t_x = diff(u_t, x); u_t_y = diff(u_t, y); u_t_z = diff(u_t, z);
v_t_x = diff(v_t, x); v_t_y = diff(v_t, y); v_t_z = diff(v_t, z);
w_t_x = diff(w_t, x); w_t_y = diff(w_t, y); w_t_z = diff(w_t, z);

% Compute grad dot pressure
p_t_x = diff(p_t, x); p_t_y = diff(p_t, y); p_t_z = diff(p_t, z);

p_t_x = simplify(p_t_x);
p_t_y = simplify(p_t_y);
p_t_z = simplify(p_t_z);

% -------------------------------------------------------------------------
% Navier-Stokes force
fx = rho0 * ( u_t + u * u_x + v * u_y + w * u_z ) + p_x - mu * (u_xx + u_yy + u_zz);

fy = rho0 * ( v_t + u * v_x + v * v_y + w * v_z ) + p_y - mu * (v_xx + v_yy + v_zz);

fz = rho0 * ( w_t + u * w_x + v * w_y + w * w_z ) + p_z - mu * (w_xx + w_yy + w_zz);
% -------------------------------------------------------------------------
% Stokes force
%fx = rho0 * u_t + p_x - mu * (u_xx + u_yy + u_zz);

%fy = rho0 * v_t + p_y - mu * (v_xx + v_yy + v_zz);

%fz = rho0 * w_t + p_z - mu * (w_xx + w_yy + w_zz);

% simplify and divide by density
fx = simplify(fx); fx = fx / rho0;
fy = simplify(fy); fy = fy / rho0;
fz = simplify(fz); fz = fz / rho0;

% traction
H = mu * [    2*u_x,   u_y+v_x,   u_z+w_x; 
            v_x+u_y,     2*v_y, v_z + w_y; 
          w_x + u_z, w_y + v_z,     2*w_z ] - p * eye(3);

H_top = H * [0;  0;  1]; H_top = simplify(H_top);
H_bot = H * [0;  0; -1]; H_bot = simplify(H_bot);
H_fro = H * [1;  0;  0]; H_fro = simplify(H_fro);
H_bac = H * [-1; 0;  0]; H_bac = simplify(H_bac);
H_rig = H * [0;  1;  0]; H_rig = simplify(H_rig);
H_lef = H * [0; -1;  0]; H_lef = simplify(H_lef);

% EOF