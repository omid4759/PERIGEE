% Plot the pressure error
clear all; close all; clc;

% Mesh size (max tet diameter)
ms = [0.0894744, 0.0662099, 0.0534411, 0.0416012, 0.0317831];

% Normalize by pipe radius
Rp = 0.3;
ms = ms / Rp;

xlim_max = 0.34;
xlim_min = 0.0925;

% Errors
err_names = { {'{\boldmath${v}$}$_h$', 'L_2'}, {'{\boldmath${v}$}$_h$', 'H_1'}, ...
              {'{\boldmath${\tau}$}$_h$', 'L_2'}, ...
              {'$p_h$', 'L_2'}, {'$p_h$', 'H_1'} };
            
err = zeros( length(err_names), length(ms) );
exa = [2.119215e+00, 2.286465e+01, 5.379971e+00, 1.953778e+01, 3.907556e+01];

% relative velo L2 err norm
err(1, :) = [1.280902e-03, 5.393694e-04, 2.717935e-04, 1.156253e-04, 4.589268e-05] / exa(1);

% relative velo H1 err norm
err(2, :) = [9.517951e-02, 5.498583e-02, 3.645852e-02, 2.290816e-02, 1.344689e-02] / exa(2);

% relative wss L2 err norm
err(3, :) = [1.861533e-02, 1.197991e-02, 9.308278e-03, 6.199857e-03, 3.827384e-03] / exa(3);

% relative pres L2 err norm
err(4, :) = [9.254529e-04, 3.690180e-04, 2.015108e-04, 1.016497e-04, 4.852524e-05] / exa(4);

% relative pres H1 err norm
err(5, :) = [5.841878e-02, 3.752442e-02, 2.898806e-02, 2.112283e-02, 1.505554e-02] / exa(5);


% Calculate rates
rates = zeros( length(err_names), length(ms) - 1 );
for i = 1 : length(err_names)
    for j = 1 : length(ms) - 1
        rates(i, j) = log( err(i,j+1) / err(i,j) ) / log( ms(j+1) / ms(j) );
    end
end

% =========== derivation ========
% rate_TH = log(y_top / y_bot) / log(x_top / x_bot);
% log(y_top / y_bot) = rate_TH * log(x_top / x_bot);
% y_top / y_bot = exp(rate_TH * log(x_top / x_bot));
% t_top = y_bot * exp(rate_TH * log(x_top / x_bot));
% y_bot = y_top / exp(rate_TH * log(x_top / x_bot))
% ================================

% Theoretical (TH) rates
rates_TH = [3.0, 2.0, 2.0, 2.5, 1.5];

fullfig();
for i = 1 : length(err_names)
    % Plot numerical results
    subaxis(2, 3, i, 'Spacing', 0.09, 'Padding', 0, 'Margin', 0.08, 'SpacingVert', 0.06);
    loglog(ms, err(i, :), 'k-o', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    hold on;
    
    % Plot triangle indicating theoretical rate
    x_top = ms(end-2);
    x_bot = 0.5*( ms(end-1) + ms(end-2) );
    y_top = err(i, end-2) - 0.5*(err(i, end-2) - err(i, end-1));
    y_bot = y_top / exp( rates_TH(i) * log( x_top / x_bot ) );
    loglog([x_bot, x_top], [y_bot, y_top], 'r', 'LineWidth', 1);
    loglog([x_bot, x_top], [y_bot, y_bot], 'r', 'LineWidth', 1);
    loglog([x_top, x_top], [y_bot, y_top], 'r', 'LineWidth', 1);
    ax = gca;
    ax.XAxis.Exponent = -1; xtickformat('%.1f');
    
    if i == 3
        text(x_bot + (x_top-x_bot)/3, y_bot - (y_top-y_bot)/2.5, '1', 'Color', 'r', 'FontSize', 12);
    elseif i == 4
        text(x_bot + (x_top-x_bot)/3, y_bot - (y_top-y_bot)/1.7, '1', 'Color', 'r', 'FontSize', 12);
    else
        text(x_bot + (x_top-x_bot)/3, y_bot - (y_top-y_bot)/2.0, '1', 'Color', 'r', 'FontSize', 12);
    end
    text(x_top + (x_top-x_bot)/3, y_bot + (y_top-y_bot)/2.5, num2str( rates_TH(i) ), 'Color', 'r', 'FontSize', 12);
    
    xlim([xlim_min, xlim_max]);
    y_range = err(i, 1) - err(i, end);
    
    if i == 1
        ylim([9.9e-6, 1.5e-3]);
    elseif i == 2
        ylim([4.0e-4, 7.0e-3]);
    elseif i == 3
        ylim([5.0e-4, 5.3e-3]); ytickformat('%.1f');
    elseif i == 4
        ylim([1.0e-6, 1.0e-4]); 
    elseif i == 5
        ylim([2.5e-4, 2.2e-3]); ytickformat('%.1f');
    end
    
    set( gca, 'Box', 'on', 'TickDir'     , 'out', ...
        'TickLength'  , [.02 .02], ...
        'XMinorTick'  , 'on'      , ...
        'YMinorTick'  , 'on'  , ...
        'YGrid'       , 'on' , ...
        'XGrid'       , 'on' , ...
        'XColor'      , [0 0 0 ], ...
        'YColor'      , [0 0 0 ], ...
        'LineWidth'   , 1 );
    axis square;
    set(gca, 'FontSize', 12,'fontWeight','bold');

    hXLabel = xlabel('$h$ / $R$', 'interpreter', 'latex');
    hYLabel = ylabel(['Relative Error of ', err_names{i}{1}, ' in $', err_names{i}{2}, '$-norm'], ...
                     'interpreter', 'latex');
    set([hXLabel, hYLabel], 'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'bold');
end

set(gcf, 'PaperOrientation','landscape');
print -dpdf conv_rate_p2.pdf -r0 -fillpage
