function compare_velo_profiles(sim_dir, z_coord, mu, rho, R, c_n, B_n, Q_n, G_n, T, n_modes, t_steps, start_step, stop_step)

colors = [     0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; ...
          0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; ...
          0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];

dt = T / t_steps;
      
omega = 2 * pi / T;                                 % base angular frequency

x = -R : 0.01 : R;  % radius vector
r_exact = abs(x);

w_fig = figure;   v_fig = figure;
w_lim = [-5, 25]; v_lim = [-7e-3, 7e-3];

sim_steps = start_step + (0 : t_steps) * (stop_step - start_step) / t_steps;

for ii = 1 : (t_steps + 1)
    

    % Assemble simulation results
    filename = [sim_dir, '/SOL_9', sprintf('%08d', sim_steps(ii)), '_velo_z7d5.txt'];
    disp(['Reading ', filename]);
    
    velo_numer = readmatrix(filename);
    
    x_numer = velo_numer(:, 1);
    y_numer = velo_numer(:, 2);
    theta_numer = atan2(y_numer, x_numer);
    r_numer = sqrt(x_numer.^2 + y_numer.^2);
    
    % Flip the radius sign for y < 0
    r_numer(y_numer < 0) = -r_numer(y_numer < 0);
    
    u_numer = velo_numer(:, 4);
    v_numer = velo_numer(:, 5);
    w_numer = velo_numer(:, 6);
    
    % Cartesian to polar transformation
    vr_numer = cos(theta_numer) .* u_numer + sin(theta_numer) .* v_numer;
    vt_numer = -sin(theta_numer) .* u_numer + cos(theta_numer) .* v_numer;
    
    t = (ii - 1) * dt;
    w_ax = subplot((t_steps + 1) / 2, 2, ii, 'Parent', w_fig); hold(w_ax, 'on');
    v_ax = subplot((t_steps + 1) / 2, 2, ii, 'Parent', v_fig); hold(v_ax, 'on');
    
    ax_all = {w_ax, v_ax};
    
    w_exact = 2 * Q_n(1) * (R^2 - r_exact.^2) / (pi * R^4);         % axial velocity

    vr_exact = zeros(1, length(r_exact));                           % radial velocity

    for k = 2 : n_modes
        
        n = k - 1;
        Omega_n  = R * sqrt(rho * n * omega / mu);   % Womersley number
        Lambda_n = 1j^1.5 * Omega_n;

        w_exact = w_exact + B_n(k) / ( rho * c_n(k) ) * ...
            ( 1 - G_n(k) * besselj(0, Lambda_n * r_exact / R) / besselj(0, Lambda_n) ) * ...
            exp(1j * n * omega * (t - z_coord / c_n(k)) );
        vr_exact = vr_exact + B_n(k) * 1j * n * omega * R / ( 2 * rho * c_n(k)^2 ) .* ...
            ( r_exact / R - G_n(k) * 2 * besselj(1, Lambda_n * r_exact / R) / (Lambda_n * besselj(0, Lambda_n)) ) * ...
            exp(1j * n * omega * (t - z_coord / c_n(k)) );
        
        plot(w_ax, real(w_exact), x / R, 'Color', colors(1, :), 'Linestyle', '-');
        plot(w_ax, w_numer, r_numer / R, 'Color', colors(2, :), 'Linestyle', 'None', 'Marker', 'o');
        plot(w_ax, [0, 0], [-1, 1], 'Color', [0.75, 0.75, 0.75] );               % grey
        
        plot(v_ax, real(vr_exact), x / R, 'Color', colors(1, :), 'Linestyle', '-' );
        plot(v_ax, vr_numer, r_numer / R, 'Color', colors(2, :), 'Linestyle', 'None', 'Marker', 'o');
        plot(v_ax, [0, 0], [-1, 1], 'Color', [0.75, 0.75, 0.75] );               % grey
        
    end
    
    w_ax.XLim  = w_lim; v_ax.XLim  = v_lim;
    
    for jj = 1 : length(ax_all)
        if ii == 1
            ylabel(ax_all{jj}, '$t$ = $0$', 'interpreter', 'latex');
        elseif ii == 2
            ylabel(ax_all{jj}, ['$t$ = $T$ / ', num2str(t_steps)], 'interpreter', 'latex');
        elseif ii == t_steps + 1
            ylabel(ax_all{jj}, '$t$ = $T$', 'interpreter', 'latex');
        else
            ylabel(ax_all{jj}, ['$t$ = ', num2str(ii-1), ' $T$ / ', num2str(t_steps)], 'interpreter', 'latex');
        end
    end
end

sgtitle(w_fig, ['Axial Velocity $w(r, z=',  num2str(z_coord), ', t)$'], 'interpreter', 'latex');
sgtitle(v_fig, ['Radial Velocity $v(r, z=', num2str(z_coord), ', t)$'], 'interpreter', 'latex');