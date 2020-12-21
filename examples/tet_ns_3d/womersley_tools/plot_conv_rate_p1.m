clear all; close all; clc;

% Mesh size (max tet diameter)
ms = [0.0212408, 0.0168113, 0.0136046, 0.0108558, 0.00926999];

% Normalize by pipe radius
Rp = 0.3;
ms = ms / Rp;

xlim_max = 0.08;
xlim_min = 0.0275;

% Errors
err_names = { {'{\boldmath${v}$}$_h$', 'L_2'}, {'{\boldmath${v}$}$_h$', 'H_1'}, ...
              {'{\boldmath${\tau}$}$_h$', 'L_2'}, ...
              {'$p_h$', 'L_2'}, {'$p_h$', 'H_1'} };
            
err = zeros( length(err_names), length(ms) );

% Exact solutions integrated over...
% exa = [2.457253e+00, 2.050088e+01, 3.086559e+00, 9.091145e-01, 5.326925e+00];  % finest p1 mesh
exa = [2.457113e+00, 2.050023e+01, 3.024062e+00, 9.091392e-01, 5.327069e+00];    % finest p2 mesh

% relative velo L2 err norm
err(1, :) = [2.699528e-03, 1.595819e-03, 1.016854e-03, 6.773337e-04, 5.009034e-04] / exa(1);

% relative velo H1 err norm
err(2, :) = [7.195206e-01, 5.686614e-01, 4.554711e-01, 3.695009e-01, 3.143580e-01] / exa(2);

% relative wss L2 err norm
err(3, :) = [1.509329e-01, 1.213661e-01, 9.932551e-02, 8.016511e-02, 6.752238e-02] / exa(3);

% relative pres L2 err norm
err(4, :) = [2.059793e-02, 1.438822e-02, 1.011933e-02, 7.320857e-03, 5.731659e-03] / exa(4);

% relative pres H1 err norm
err(5, :) = [2.421704e+00, 2.127668e+00, 1.883097e+00, 1.683841e+00, 1.548086e+00] / exa(5);

% Calculate rates
rates = zeros( length(err_names), length(ms) - 1 );
for i = 1 : length(err_names)
    for j = 1 : length(ms) - 1
        rates(i, j) = log( err(i,j+1) / err(i,j) ) / log( ms(j+1) / ms(j) );
    end
end

% Theoretical (TH) rates
rates_TH = [2.0, 1.0, 1.0, 1.5, 0.5];

fullfig();
for i = 1 : length(err_names)

    subaxis(2, 3, i, 'Spacing', 0.09, 'Padding', 0, 'Margin', 0.08, 'SpacingVert', 0.06);
    loglog(ms, err(i, :), 'k-o', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    hold on;
    
% =========== Plot triangle indicating theoretical rate ========
% rate_TH = log(y_top / y_bot) / log(x_top / x_bot);
% log(y_top / y_bot) = rate_TH * log(x_top / x_bot);
% y_top / y_bot = exp(rate_TH * log(x_top / x_bot));
% t_top = y_bot * exp(rate_TH * log(x_top / x_bot));
% y_bot = y_top / exp(rate_TH * log(x_top / x_bot))
% ==============================================================
%     x_top = ms(end-2);
%     x_bot = 0.5*( ms(end-1) + ms(end-2) );
%     y_top = err(i, end-2) - 0.5*(err(i, end-2) - err(i, end-1));
%     y_bot = y_top / exp( rates_TH(i) * log( x_top / x_bot ) );
%     loglog([x_bot, x_top], [y_bot, y_top], 'r', 'LineWidth', 1);
%     loglog([x_bot, x_top], [y_bot, y_bot], 'r', 'LineWidth', 1);
%     loglog([x_top, x_top], [y_bot, y_top], 'r', 'LineWidth', 1);
%
%     if i == 1
%         text(x_bot + (x_top-x_bot)/3, y_bot - (y_top-y_bot)/2.0, '1', 'Color', 'r', 'FontSize', 12);
%     else
%         text(x_bot + (x_top-x_bot)/3, y_bot - (y_top-y_bot)/1.5, '1', 'Color', 'r', 'FontSize', 12);
%     end
%     text(x_top + (x_top-x_bot)/3, y_bot + (y_top-y_bot)/2.5, num2str( rates_TH(i) ), 'Color', 'r', 'FontSize', 12);

% ========== Annotate rates ========
% 2 * log(x) = log(ms(j+1)) + log(ms(j)) = log( ms(j+1)*ms(j) )
% log(x) = log( ms(j+1)*ms(j) ) / 2
% x = exp(log( ms(j+1)*ms(j) ) / 2)
% ===================================
    for j = 1 : length(ms) - 1
        text( exp( log(ms(j+1)*ms(j)) / 2.03), exp( log(err(i,j+1)*err(i,j)) / 1.98 ), ...
               num2str( rates(i, j), '%.2f'), 'Color', [0, 0, 0.545], 'FontSize', 12);
    end
    
    ax = gca;
    ax.XAxis.Exponent = -2;
    
    xlim([xlim_min, xlim_max]);
    y_range = err(i, 1) - err(i, end);
    
    if i == 1
        ylim([1.70e-4, 1.30e-3]); ax.YAxis.Exponent = -3;
    elseif i == 2
        ylim([1.40e-2, 3.90e-2]); ax.YAxis.Exponent = -2;
    elseif i == 3
        ylim([2.02e-2, 5.63e-2]); ax.YAxis.Exponent = -2;
    elseif i == 4
        ylim([5.50e-3, 2.60e-2]); ax.YAxis.Exponent = -2;
    elseif i == 5
        ylim([2.78e-1, 4.80e-1]); ax.YAxis.Exponent = -1; % ytickformat('%.1f');
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
    set(gca, 'FontSize', 12, 'fontWeight','bold');

    hXLabel = xlabel('$h$ / $R$', 'interpreter', 'latex');
    hYLabel = ylabel(['Relative Error of ', err_names{i}{1}, ' in $', err_names{i}{2}, '$-norm'], ...
                     'interpreter', 'latex');
    set([hXLabel, hYLabel], 'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'bold');
end

set(gcf, 'PaperOrientation','landscape');
print -dpdf conv_rate_p1.pdf -r0 -fillpage
