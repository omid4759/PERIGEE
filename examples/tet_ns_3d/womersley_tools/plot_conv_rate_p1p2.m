% 12/20/2020 - Only plot 2x2 grid with velo L2, wss L2, pres L2, pres H1,
% since velo H1 is very similar to wss L2

clear all; close all; clc;

% Mesh size (max tet diameter)
ms = [0.0212408, 0.0168113, 0.0136046, 0.0108558, 0.00926999;   % p1
      0.0412092, 0.0335024, 0.0269549, 0.0215289, 0.0185061  ]; % p2

% Normalize by pipe radius
Rp = 0.3;
ms = ms / Rp;

xlim_max = 0.155;
xlim_min = 0.0275;

% Errors
% err_names = { {'{\boldmath${v}$}$_h$', 'L_2'}, {'{\boldmath${v}$}$_h$', 'H_1'}, ...
%               {'{\boldmath${\tau}$}$_h$', 'L_2'}, ...
%               {'$p_h$', 'L_2'}, {'$p_h$', 'H_1'} };
err_names = { {'{\boldmath${v}$}$_h$', 'L_2'}, {'{\boldmath${\tau}$}$_h$', 'L_2'}, ...
              {'$p_h$', 'L_2'}, {'$p_h$', 'H_1'} };
          
err = zeros( length(err_names), length(ms), 2 );

% Exact solutions integrated over the finest P2 mesh
% exa = [2.457113e+00, 2.050023e+01, 3.024062e+00, 9.091392e-01, 5.327069e+00];
exa = [2.457113e+00, 3.024062e+00, 9.091392e-01, 5.327069e+00];

% relative velo L2 err norm
err(1, :, :) = [2.699528e-03, 1.595819e-03, 1.016854e-03, 6.773337e-04, 5.009034e-04; 
                1.183015e-04, 6.582388e-05, 3.231378e-05, 1.748338e-05, 1.064649e-05 ]' / exa(1);

% % relative velo H1 err norm
% err(2, :, :) = [7.195206e-01, 5.686614e-01, 4.554711e-01, 3.695009e-01, 3.143580e-01;
%                 2.556999e-02, 1.606325e-02, 1.022238e-02, 6.623977e-03, 4.917036e-03 ]' / exa(2);

% relative wss L2 err norm
err(2, :, :) = [1.509329e-01, 1.213661e-01, 9.932551e-02, 8.016511e-02, 6.752238e-02;
                6.833765e-03, 4.340224e-03, 2.773377e-03, 1.784367e-03, 1.327418e-03 ]' / exa(2);


% relative pres L2 err norm
err(3, :, :) = [2.059793e-02, 1.438822e-02, 1.011933e-02, 7.320857e-03, 5.731659e-03;
                5.006683e-04, 2.284029e-04, 1.079442e-04, 5.261179e-05, 3.309989e-05 ]' / exa(3);

% relative pres H1 err norm
err(4, :, :) = [2.421704e+00, 2.127668e+00, 1.883097e+00, 1.683841e+00, 1.548086e+00;
                3.757109e-02, 2.494294e-02, 1.705044e-02, 1.203480e-02, 9.495378e-03 ]' / exa(4);

% Calculate rates
rates = zeros( length(err_names), length(ms)-1, 2 );
for k = 1 : 2
    for i = 1 : length(err_names)
        for j = 1 : length(ms) - 1
            rates(i,j,k) = log( err(i,j+1,k) / err(i,j,k) ) / log( ms(k,j+1) / ms(k,j) );
        end
    end
end

% Theoretical (TH) rates
% rates_TH = [2.0, 1.0, 1.0, 1.5, 0.5;
%             3.0, 2.0, 2.0, 3.0, 1.5 ];
rates_TH = [2.0, 1.0, 1.5, 0.5;
            3.0, 2.0, 3.0, 1.5 ];

color  = [0.918, 0.235, 0.325;
              0,     0, 0.545 ];
marker = ['o', 's'];


% Precise control over position of text annotations
x_annot_loc(:, :, 1) = [2.05, 2.04, 2.02, 2.01;
                        2.06, 2.05, 2.03, 1.99;
                        2.01, 2.01, 2.00, 1.98;
                        2.04, 2.02, 2.00, 1.96];
x_annot_loc(:, :, 2) = [1.79, 1.82, 1.84, 1.85;
                        1.81, 1.82, 1.84, 1.82;
                        1.81, 1.81, 1.83, 1.81;
                        1.80, 1.81, 1.83, 1.82];
y_annot_loc(:, :, 1) = [1.94, 1.95, 1.94, 1.93;
                        1.87, 1.87, 1.86, 1.85;
                        1.84, 1.84, 1.84, 1.81;
                        1.42, 1.44, 1.44, 1.43];
y_annot_loc(:, :, 2) = [2.04, 2.04, 2.05, 2.04;
                        2.06, 2.06, 2.06, 2.04;
                        2.07, 2.07, 2.07, 2.05;
                        2.11, 2.11, 2.10, 2.07];
figure;
for i = 1 : length(err_names)

    subaxis(2, 2, i, 'Spacing', 0.15, 'Padding', 0, 'MarginLeft', 0.1, ...
            'MarginRight', 0.02, 'SpacingVert', 0.06);
    for k = 1 : 2
        loglog(ms(k, :), err(i, :, k), '-', 'LineWidth', 1, 'Marker', marker(k), 'MarkerSize', 4, ...
               'Color', color(k,:), 'MarkerFaceColor', color(k,:)); hold on;
        
        % Annotate all rates
        for j = 1 : length(ms) - 1
            text( exp( log(ms(k,j+1)*ms(k,j)) / x_annot_loc(i,j,k)), ...
                  exp( log(err(i,j+1,k)*err(i,j,k)) / y_annot_loc(i,j,k) ), ...
                  num2str( rates(i,j,k), '%.2f'), 'Color', color(k,:), 'FontSize', 12);
        end
        
    end
    
    ax = gca;
    ax.XAxis.Exponent = -1;
    
    xlim([xlim_min, xlim_max]);
    
    if i == 1
        ylim([3.00e-6, 1.60e-3]);
    elseif i == 2
        ylim([3.3e-4, 6.90e-2]); 
    elseif i == 3
        ylim([2.40e-5, 3.70e-2]); 
    elseif i == 4
        ylim([1.2e-3, 7.00e-1]);
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

lg = legend('P1', 'P2', 'Orientation', 'horizontal', 'Box', 'off');
set(lg, 'Position', [0.4, -0.08, 0.2, 0.2], 'Units', 'normalized');

set(gcf, 'WindowState','fullscreen');
% set(gcf, 'PaperOrientation','landscape');
print -dpdf conv_rate_p1p2.pdf -r0 -fillpage