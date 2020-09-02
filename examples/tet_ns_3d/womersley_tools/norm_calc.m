% calculate the norm of quantities
clear all; clc;

t = 1.0;
mu = 0.1;
rho0 = 1.0;
a = pi / 4;
d = pi / 2;

syms x y z;

% Given analytic manufactured solution
u = -a * (exp(a*x) * sin(a*y + d*z) + exp(a*z) * cos(a*x + d*y)) * ...
    exp(-mu/rho0 * d^2 * t);

v = -a * (exp(a*y) * sin(a*z + d*x) + exp(a*x) * cos(a*y + d*z)) * ...
    exp(-mu/rho0 * d^2 * t);

w = -a * (exp(a*z) * sin(a*x + d*y) + exp(a*y) * cos(a*z + d*x)) * ...
    exp(-mu/rho0 * d^2 * t);

p = -a^2 / 2 * (exp(2*a*x) + exp(2*a*y) + exp(2*a*z) + ...
    2 * sin(a*x + d*y) * cos(a*z + d*x) * exp(a*(y + z)) + ...
    2 * sin(a*y + d*z) * cos(a*x + d*y) * exp(a*(z + x)) + ...
    2 * sin(a*z + d*x) * cos(a*y + d*z) * exp(a*(x + y))) * ...
    exp(-2 * mu/rho0 * d^2 * t);

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

% -------------------------------------------------------------------------
v_l2 = double( sqrt( int( int( int( u*u+v*v+w*w, x, -1, 1), y, -1, 1), z, -1, 1) ) );

v_h1 = double( sqrt( int( int( int( u*u + v*v + w*w ...
    + u_x*u_x + v_x*v_x + w_x*w_x ...
    + u_y*u_y + v_y*v_y + w_y*w_y ...
    + u_z*u_z + v_z*v_z + w_z*w_z, x, -1, 1), y, -1, 1), z, -1, 1) ) );

p_l2 = double( sqrt( int( int( int( p*p, x, -1, 1), y, -1, 1), z, -1, 1) ) );

% eof