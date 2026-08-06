%close all;

% -------------------------------------------------------------------------
% Controller / estimator settings
% -------------------------------------------------------------------------
Q = sys{11}.C' * sys{11}.C;        % State penalty regulator
Qi = 5e-7;
Qi = blkdiag(Q, Qi);               % Integral action penalty

W = 1;                             % Process covariance Kalman filter
V = 3e-4;                          % Noise covariance Kalman filter

R = [5e-6 75e-7 1e-5 3e-5 5e-5];
% Frequency vector
w = logspace(-5, 3, 30000);

% -------------------------------------------------------------------------
% Create figure and main axes
% -------------------------------------------------------------------------
figure
ax = axes;
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');

% Set colour and line-style order before plotting
colororder(ax, [ ...
    0/255   178/255 114/255; ...
    230/255 159/255   0/255; ...
    0/255     0/255   0/255; ...
    0/255   121/255 167/255; ...
    240/255 228/255  66/255; ...
    0/255   158/255 115/255; ...
    86/255  180/255 233/255; ...
    213/255  94/255   0/255 ...
]);

%linestyles = ["-.","--","-",":","-."];
%linestyleorder(ax, linestyles);

labels = { ...
    '$R=5\times 10^{-6}$', ...
    '$R=7.5\times 10^{-6}$', ...
    '$R=1\times 10^{-5}$', ...
    '$R=3\times 10^{-5}$', ...
    '$R=5\times 10^{-5}$'};

% Compute Strouhal number from frequency output
[~,~,wout] = bode(sys{idx_list(1)}, w);
St = wout(:) / (U_inf * 2*pi);

% Store data for inset
magdb_all = cell(1, numel(R));

% -------------------------------------------------------------------------
% Plot R sweep on main axes
% -------------------------------------------------------------------------
for i = 1:numel(R)

    [~, ~, L_sweep{i}, ~, R_sweep_H{i}, R_sweep_H2{i}, R_sweep_H3{i}, ~, ~] = ...
        Sys_response(sys{11}.A, sys{11}.B, CLQR{11}, CKF{11}, ...
        Ts, Q, Qi, R(i), W, V, r);
    [Gm{i}, Pm{i}, ~, ~] = margin(L_sweep{i});
    [mag, ~, ~] = bode(R_sweep_H3{i}, w);
    mag = squeeze(mag);
    magdb = mag2db(mag);

    magdb_all{i} = magdb;

    semilogx(ax, St, magdb, ...
        'DisplayName', labels{i}, ...
        'LineWidth', 1.5);
end

% -------------------------------------------------------------------------
% Main axes labels and limits
% -------------------------------------------------------------------------
xlabel(ax, 'Normalised Frequency, $\mathrm{St}$', ...
    'Interpreter', 'latex', ...
    'FontSize', 16);

ylabel(ax, '$Magnitude, |\tilde{G}(e^{j\varphi})|$ (dB)', ...
    'Interpreter', 'latex', ...
    'FontSize', 16);

xlim(ax, [5e-4 3]);

legend(ax, 'show', ...
    'FontSize', 12, ...
    'Location', 'best', ...
    'Interpreter', 'latex');

% -------------------------------------------------------------------------
% Zoom region for inset
% -------------------------------------------------------------------------
x_zoom = [0.01 2];
y_zoom = [-3 1.5];

% Add zoom rectangle to main axes
% Use line segments instead of rectangle because this behaves better on log axes.
plot(ax, [x_zoom(1) x_zoom(2)], [y_zoom(1) y_zoom(1)], ...
    'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

plot(ax, [x_zoom(1) x_zoom(2)], [y_zoom(2) y_zoom(2)], ...
    'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

plot(ax, [x_zoom(1) x_zoom(1)], [y_zoom(1) y_zoom(2)], ...
    'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

plot(ax, [x_zoom(2) x_zoom(2)], [y_zoom(1) y_zoom(2)], ...
    'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

% Re-apply main limits after drawing rectangle
xlim(ax, [5e-4 3]);

% -------------------------------------------------------------------------
% Create inset axes
% -------------------------------------------------------------------------
ax2 = axes('Position', [0.3 0.2 0.32 0.32]);
hold(ax2, 'on');
grid(ax2, 'on');
box(ax2, 'on');

colororder(ax2, colororder(ax));
linestyleorder(ax2, linestyles);

% Logical index for zoomed region
Lv = (St >= x_zoom(1)) & (St <= x_zoom(2));

% Plot zoomed data
for i = 1:numel(R)
    semilogx(ax2, St(Lv), magdb_all{i}(Lv), ...
        'LineWidth', 1.5);
end

xlim(ax2, x_zoom);
ylim(ax2, y_zoom);

% Inset tick formatting
set(ax2, 'XScale', 'log');
ax2.FontSize = 12;
ax2.XTick = [];
ax2.YTick = [-4 0 2];
ax2.YTickLabel = {'-4', '0', '2'};

% -------------------------------------------------------------------------
% Figure formatting
% -------------------------------------------------------------------------
set(gcf, 'Position', [50 50 1300 500]);

set(ax, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 16, ...
    'TickDir', 'in', ...
    'Box', 'off', ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on');

set(ax, 'XScale', 'log');
xlim(ax, [5e-4 3]);

set(ax2, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 12, ...
    'TickDir', 'in', ...
    'Box', 'on', ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on');

% Main x-axis tick formatting
set(ax, 'XTick', [1e-4 1e-3 1e-2 1e-1 1]);
set(ax, 'XTickLabel', { ...
    '$10^{-4}$', ...
    '$10^{-3}$', ...
    '$10^{-2}$', ...
    '$10^{-1}$', ...
    '$1$'});

ax.TickLabelInterpreter = 'latex';

% -------------------------------------------------------------------------
% Export
% -------------------------------------------------------------------------
if print_flag == 1
    exportgraphics(gcf, 'R_sweep_bode_inset.pdf', 'ContentType', 'vector');
    close all;
end
%%
%close all;

% -------------------------------------------------------------------------
% Controller / estimator settings
% -------------------------------------------------------------------------
Q = sys{11}.C' * sys{11}.C;        % State penalty regulator
Qi = 5e-7;
Qi = blkdiag(Q, Qi);               % Integral action penalty

q=[1e-9 1e-8 1e-7 1e-6 1e-2 1e-0];

V = 1e-4;

% Frequency vector
w = logspace(-5, 3, 30000);

% -------------------------------------------------------------------------
% Create figure and main axes
% -------------------------------------------------------------------------
figure
ax = axes;
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');

% Set colour and line-style order before plotting
colororder(ax, [ ...
    0/255   178/255 114/255; ...
    230/255 159/255   0/255; ...
    0/255     0/255   0/255; ...
    0/255   121/255 167/255; ...
    240/255 228/255  66/255; ...
    0/255   158/255 115/255; ...
    86/255  180/255 233/255; ...
    213/255  94/255   0/255 ...
]);

%linestyles = ["-.","--","-",":","-."];
%linestyleorder(ax, linestyles);

labels = { ...
    '$W/V=1\times 10^{3}$', ...
    '$W/V=1\times 10^{4}$', ...
    '$W/V=3.3\times 10^{4}$', ...
    '$W/V=1\times 10^{5}$', ...
    '$W/V=1\times 10^{6}$',...
    '$W/V=1\times 10^{7}$',...
    'LQR'};

% Compute Strouhal number from frequency output
[~,~,wout] = bode(sys{idx_list(1)}, w);
St = wout(:) / (U_inf * 2*pi);

% Store data for inset
magdb_all = cell(1, numel(V));

for i = 1:numel(q)
    W = eye(11)+q(i)*sys{11}.B*sys{11}.B';
    Wi = eye(12)+blkdiag(q(i)*sys{11}.B*sys{11}.B',0);
    [~, ~, L_sweep{i}, L_lqri{i}, R_sweep_H{i}, R_sweep_H2{i}, R_sweep_H3{i}, ~, ~] = ...
        Sys_response(sys{11}.A, sys{11}.B, CLQR{11}, CKF{11}, ...
        Ts, Q, Qi, R(i), W, Wi, V, r);
    [Gm{i}, Pm{i}, ~, ~] = margin(L_sweep{i});
    [mag, ~, ~] = bode(R_sweep_H3{i}, w);
    mag = squeeze(mag);
    magdb = mag2db(mag);

    magdb_all{i} = magdb;


    semilogx(ax, St, magdb, ...
        'DisplayName', labels{i}, ...
        'LineWidth', 1.5);
end
i=i+1;
H3_lqri = minreal(1/((L_lqri{1})+1));
[mag, ~, ~] = bode(H3_lqri, w);
mag = squeeze(mag);
magdb = mag2db(mag);
semilogx(ax, St, magdb, 'DisplayName', labels{i}, 'LineWidth', 1.5);
% -------------------------------------------------------------------------
% Main axes labels and limits
% -------------------------------------------------------------------------
xlabel(ax, 'Normalised Frequency, $\mathrm{St}$', ...
    'Interpreter', 'latex', ...
    'FontSize', 16);

ylabel(ax, '$Magnitude, |\tilde{G}(e^{j\varphi})|$ (dB)', ...
    'Interpreter', 'latex', ...
    'FontSize', 16);

xlim(ax, [5e-4 3]);

legend(ax, 'show', ...
    'FontSize', 12, ...
    'Location', 'best', ...
    'Interpreter', 'latex');
%{
% -------------------------------------------------------------------------
% Zoom region for inset
% -------------------------------------------------------------------------
x_zoom = [0.01 2];
y_zoom = [-3 1.5];

% Add zoom rectangle to main axes
% Use line segments instead of rectangle because this behaves better on log axes.
plot(ax, [x_zoom(1) x_zoom(2)], [y_zoom(1) y_zoom(1)], ...
    'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

plot(ax, [x_zoom(1) x_zoom(2)], [y_zoom(2) y_zoom(2)], ...
    'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

plot(ax, [x_zoom(1) x_zoom(1)], [y_zoom(1) y_zoom(2)], ...
    'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

plot(ax, [x_zoom(2) x_zoom(2)], [y_zoom(1) y_zoom(2)], ...
    'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

% Re-apply main limits after drawing rectangle
xlim(ax, [5e-4 3]);

% -------------------------------------------------------------------------
% Create inset axes
% -------------------------------------------------------------------------
ax2 = axes('Position', [0.3 0.2 0.32 0.32]);
hold(ax2, 'on');
grid(ax2, 'on');
box(ax2, 'on');

colororder(ax2, colororder(ax));
linestyleorder(ax2, linestyles);

% Logical index for zoomed region
Lv = (St >= x_zoom(1)) & (St <= x_zoom(2));

% Plot zoomed data
for i = 1:numel(R)
    semilogx(ax2, St(Lv), magdb_all{i}(Lv), ...
        'LineWidth', 1.5);
end

xlim(ax2, x_zoom);
ylim(ax2, y_zoom);

% Inset tick formatting
set(ax2, 'XScale', 'log');
ax2.FontSize = 12;
ax2.XTick = [];
ax2.YTick = [-4 0 2];
ax2.YTickLabel = {'-4', '0', '2'};
%}
% -------------------------------------------------------------------------
% Figure formatting
% -------------------------------------------------------------------------
set(gcf, 'Position', [50 50 1300 500]);

set(ax, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 16, ...
    'TickDir', 'in', ...
    'Box', 'off', ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on');

set(ax, 'XScale', 'log');
xlim(ax, [5e-4 3]);
%{
set(ax2, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 12, ...
    'TickDir', 'in', ...
    'Box', 'on', ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on');

% Main x-axis tick formatting
set(ax, 'XTick', [1e-4 1e-3 1e-2 1e-1 1]);
set(ax, 'XTickLabel', { ...
    '$10^{-4}$', ...
    '$10^{-3}$', ...
    '$10^{-2}$', ...
    '$10^{-1}$', ...
    '$1$'});

ax.TickLabelInterpreter = 'latex';

% -------------------------------------------------------------------------
% Export
% -------------------------------------------------------------------------
if print_flag == 1
    exportgraphics(gcf, 'R_sweep_bode_inset.pdf', 'ContentType', 'vector');
    close all;
end
%}
%%

Q = sys{11}.C' * sys{11}.C;        % State penalty regulator
Qi = 5e-7;
Qi = blkdiag(Q, Qi);               % Integral action penalty

W = 1;                             % Process covariance Kalman filter             
R = logspace(-8,-1,30);
V = logspace(-8,-1,30);

for i = 1:numel(R)
    for j = 1:numel(V)

        [~, ~, L_sweep{i,j}, ~, ~, ~, ~, ~, ~] = ...
        Sys_response(sys{11}.A, sys{11}.B, CLQR{11}, CKF{11}, ...
        Ts, Q, Qi, R(i), W, V(j), r);
        [Gm{i,j}, Pm{i,j}, ~, ~] = margin(L_sweep{i,j});
    end
end
%%

Z = 20*log10(cell2mat(Gm));
Z2 =cell2mat(Pm);
[X,Y] = meshgrid(R,V);
figure;
surf(X,Y,Z,20)
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
axis tight; % Adjust the axis limits 
xlabel('R')
ylabel('V')
cb = colorbar(); 
ylabel(cb,'Magnitude (dB)','FontSize',16,'Rotation',270)
figure
surf(X,Y,Z2,20)
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
axis tight; % Adjust the axis limits 
xlabel('R')
ylabel('V')
cb = colorbar(); 
ylabel(cb,'Phase (deg)','FontSize',16,'Rotation',270)