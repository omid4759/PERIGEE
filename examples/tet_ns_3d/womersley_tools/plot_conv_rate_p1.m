% Plot the pressure error
clear all; close all; clc;

% Mesh size (max tet diameter)
ms = [0.0412092, 0.0335024, 0.0271552, 0.0212408, 0.0160842];

% Normalize by pipe radius
Rp = 0.3;
ms = ms / Rp;

xlim_max = 0.155;
xlim_min = 0.0475;

% Errors
err_names = { {'{\boldmath${v}$}$_h$', 'L_2'}, {'{\boldmath${v}$}$_h$', 'H_1'}, ...
              {'{\boldmath${\tau}$}$_h$', 'L_2'}, ...
              {'$p_h$', 'L_2'}, {'$p_h$', 'H_1'} };
            
err = zeros( length(err_names), length(ms) );
exa = [2.118155e+00, 2.284145e+01, 5.190399e+00, 1.953617e+01, 3.907234e+01];

% relative velo L2 err norm
err(1, :) = [1.501688e-02, 1.009462e-02, 6.836738e-03, 4.338830e-03, 2.545145e-03] / exa(1);

% relative velo H1 err norm
err(2, :) = [6.709449e-01, 5.366545e-01, 4.280984e-01, 3.383041e-01, 2.567001e-01] / exa(2);

% relative wss L2 err norm
err(3, :) = [9.242385e-02, 7.311273e-02, 5.790696e-02, 4.363671e-02, 3.313928e-02] / exa(3);

% relative pres L2 err norm
err(4, :) = [1.069144e-02, 7.614018e-03, 5.454618e-03, 3.769337e-03, 2.497782e-03] / exa(4);

% relative pres H1 err norm
err(5, :) = [6.250368e-01, 5.548113e-01, 4.949569e-01, 4.371132e-01, 3.801927e-01] / exa(5);


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
rates_TH = [2.0, 1.0, 1.0, 1.5, 0.5];

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
    ax.XAxis.Exponent = -1;
    
    if i == 1
        text(x_bot + (x_top-x_bot)/3, y_bot - (y_top-y_bot)/2.0, '1', 'Color', 'r', 'FontSize', 12);
    else
        text(x_bot + (x_top-x_bot)/3, y_bot - (y_top-y_bot)/1.5, '1', 'Color', 'r', 'FontSize', 12);
    end
    text(x_top + (x_top-x_bot)/3, y_bot + (y_top-y_bot)/2.5, num2str( rates_TH(i) ), 'Color', 'r', 'FontSize', 12);
    
    xlim([xlim_min, xlim_max]);
    y_range = err(i, 1) - err(i, end);
    
    if i == 1
        ylim([9.0e-4, 1.0e-2]);
    elseif i == 2
        ylim([9.0e-3, 3.7e-2]); ax.YAxis.Exponent = -2; ytickformat('%.1f');
    elseif i == 3
        ylim([5.0e-3, 2.3e-2]); ax.YAxis.Exponent = -2; ytickformat('%.1f');
    elseif i == 4
        ylim([9.0e-5, 8.0e-4]); ytickformat('%.1f');
    elseif i == 5
        ylim([8.5e-3, 1.8e-2]); ax.YAxis.Exponent = -2; ytickformat('%.1f');
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
print -dpdf conv_rate_p1.pdf -r0 -fillpage
