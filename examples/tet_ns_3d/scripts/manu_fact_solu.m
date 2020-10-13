% manufactured solution for incompressible Navier-Stokes 3D
clear all; clc;

syms x y z t nu R;

u = 0;

v = 0;

% axial velo corresponding to w = r^4 - 0.09 * r^2
w = x*x*x*x + y*y*y*y + 2*x*x*y*y - 0.09*x*x - 0.09*y*y;
p = -30*z;

% compute first order derivatives
u_x = diff(u,x); u_y = diff(u,y); u_z = diff(u,z);
v_x = diff(v,x); v_y = diff(v,y); v_z = diff(v,z);
w_x = diff(w,x); w_y = diff(w,y); w_z = diff(w,z);

u_t = diff(u,t); v_t = diff(v,t); w_t = diff(w,t);

u_xx = diff(u_x, x); u_yy = diff(u_y, y); u_zz = diff(u_z, z);
v_xx = diff(v_x, x); v_yy = diff(v_y, y); v_zz = diff(v_z, z);
w_xx = diff(w_x, x); w_yy = diff(w_y, y); w_zz = diff(w_z, z);

p_x = diff(p, x); p_y = diff(p, y); p_z = diff(p, z);

div_velo = u_x + v_y + w_z; div_velo = simplify(div_velo);

% Remember to divide by rho in perigee code
fx = u_t + u * u_x + v * u_y + w * u_z + p_x - nu * (u_xx + u_yy + u_zz);

fy = v_t + u * v_x + v * v_y + w * v_z + p_y - nu * (v_xx + v_yy + v_zz);

fz = w_t + u * w_x + v * w_y + w * w_z + p_z - nu * (w_xx + w_yy + w_zz);

% traction
H = nu* [2*u_x, u_y+v_x, u_z+w_x; 
  v_x+u_y, 2*v_y, v_z + w_y; 
  w_x + u_z, w_y + v_z, 2*w_z] - p *eye(3);

H_top = H * [0;0;1]; H_top = simplify(H_top);
H_bot = H * [0;0;-1]; H_bot = simplify(H_bot);
H_fro = H * [1;0;0]; H_fro = simplify(H_fro);
H_bac = H * [-1;0;0]; H_back = simplify(H_bac);
H_rig = H * [0;1;0]; H_rig = simplify(H_rig);
H_lef = H * [0;-1;0]; H_lef = simplify(H_lef);



% EOF