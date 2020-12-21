clear all; close all; clc;

% Mesh size (max tet diameter)
ms = [0.0412092, 0.0335024, 0.0269549, 0.0215289, 0.0185061];

% Normalize by pipe radius
Rp = 0.3;
ms = ms / Rp;

xlim_max = 0.154;
xlim_min = 0.056;

% Errors
err_names = { {'{\boldmath${v}$}$_h$', 'L_2'}, {'{\boldmath${v}$}$_h$', 'H_1'}, ...
              {'{\boldmath${\tau}$}$_h$', 'L_2'}, ...
              {'$p_h$', 'L_2'}, {'$p_h$', 'H_1'} };
            
err = zeros( length(err_names), length(ms) );
exa = [2.457113e+00, 2.050023e+01, 3.024062e+00, 9.091392e-01, 5.327069e+00];

% relative velo L2 err norm
err(1, :) = [1.183015e-04, 6.582388e-05, 3.231378e-05, 1.748338e-05, 1.064649e-05] / exa(1);

% relative velo H1 err norm
err(2, :) = [2.556999e-02, 1.606325e-02, 1.022238e-02, 6.623977e-03, 4.917036e-03] / exa(2);

% relative wss L2 err norm
err(3, :) = [6.833765e-03, 4.340224e-03, 2.773377e-03, 1.784367e-03, 1.327418e-03] / exa(3);

% relative pres L2 err norm
err(4, :) = [5.006683e-04, 2.284029e-04, 1.079442e-04, 5.261179e-05, 3.309989e-05] / exa(4);

% relative pres H1 err norm
err(5, :) = [3.757109e-02, 2.494294e-02, 1.705044e-02, 1.203480e-02, 9.495378e-03] / exa(5);


% Calculate rates
rates = zeros( length(err_names), length(ms) - 1 );
for i = 1 : length(err_names)
    for j = 1 : length(ms) - 1
        rates(i, j) = log( err(i,j+1) / err(i,j) ) / log( ms(j+1) / ms(j) );
    end
end

% Theoretical (TH) rates
rates_TH = [3.0, 2.0, 2.0, 2.5, 1.5];

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
%     if i == 3
%         text(x_bot + (x_top-x_bot)/3, y_bot - (y_top-y_bot)/2.5, '1', 'Color', 'r', 'FontSize', 12);
%     elseif i == 4
%         text(x_bot + (x_top-x_bot)/3, y_bot - (y_top-y_bot)/1.7, '1', 'Color', 'r', 'FontSize', 12);
%     else
%         text(x_bot + (x_top-x_bot)/3, y_bot - (y_top-y_bot)/2.0, '1', 'Color', 'r', 'FontSize', 12);
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
    ax.XAxis.Exponent = -1;
    
    xlim([xlim_min, xlim_max]);
    y_range = err(i, 1) - err(i, end);
    
    if i == 1
        ylim([3.00e-6, 6.8e-5]);
    elseif i == 2
        ylim([1.80e-4, 1.60e-3]); ax.YAxis.Exponent = -3;
    elseif i == 3
        ylim([3.55e-4, 2.85e-3]); 
    elseif i == 4
        ylim([2.50e-5, 8.0e-4]); 
    elseif i == 5
        ylim([1.45e-3, 8.5e-3]);
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
