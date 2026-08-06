R_sweep = [5e-6 7.5e-6 1e-5 3e-5 5e-5];
for i = 1:5
    Q=(CLQR{11}'*CLQR{11}); %State penalty Regulator
    Qi = 5e-6;
    Qi = blkdiag(Q,Qi); %Integral action penalty
    R = R_sweep(i); %Actuation penalty Regulator
    W = 1; % process covariance Kalman filter
    V = 3e-4;% noise covariance Kalman filter
    [~, ~, ~, ~, ~, ~, H3_R_Sweep{i}, ~, ~] = Sys_response(sys{11}.A,sys{11}.B,CLQR{11},CKF{11},Ts,Q,Qi,R,W,V,r);
end

% Frequency vector
w = logspace(-5, 0, 30000);

% Create figure and main axes
figure
ax = axes;
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');


colororder([ ...
    0/255 0/255 0/255;    ...
    230/255 20/255   0/255;    ...
    0/255   200/255   0/255;    ...
    20/255 0 200/255;    ...
    250/255 0/255 250/255;     ...
    0/255 150/255 150/255;     ...
    ]);


[~, ~, wout] = bode(H3_R_Sweep{1}, w);
St = wout(:) / (U_inf * 2*pi);


for k = 1:numel(R_sweep)
    i = R_sweep(k);

    [mag, ~, ~] = bode(H3_R_Sweep{k}, w);
    mag = squeeze(mag);
    magdb = mag2db(mag);

    semilogx(ax, St, magdb,'LineWidth',1.5);
end


xlabel(ax, 'Strouhal number, $\mathrm{St}$', ...
    'Interpreter', 'latex', ...
    'FontSize', 16);

ylabel(ax, 'Sensitivity, $|\tilde{S}(e^{j\varphi})|$ (dB)', ...
    'Interpreter', 'latex', ...
    'FontSize', 16);
legend('R=$5 \times 10^{-6}$','R=$7.5 \times 10^{-6}$','R=$1 \times 10^{-5}$','R=$3 \times 10^{-5}$','R=$5 \times 10^{-5}$', ...
    Location='northwest')



x_zoom = [0.01 1.5];
y_zoom = [-5 5];


rectangle(ax, ...
    'Position', [x_zoom(1), y_zoom(1), diff(x_zoom), diff(y_zoom)], ...
    'LineStyle', '--', ...
    'EdgeColor', 'k', ...
    'LineWidth', 1);

xlim(ax, [1e-4 2]);
ylim(ax, [-45 8]);


ax2 = axes('Position', [0.4 0.22 0.25 0.3]);
hold(ax2, 'on');
grid(ax2, 'on');
box(ax2, 'on');

colororder(ax2, colororder(ax));

Lv = (St >= x_zoom(1)) & (St <= x_zoom(2));

% Plot zoomed data
for k = 1:numel(R_sweep)
    i = R_sweep(k);

    [mag, ~, ~] = bode(H3_R_Sweep{k}, w);
    mag = squeeze(mag);
    magdb = mag2db(mag);

    semilogx(ax2, St(Lv), magdb(Lv),'LineWidth',1.5);
end


xlim(ax2, x_zoom);
ylim(ax2, y_zoom);

ax2.XTick = [];
ax2.YTick = [-3 0];
ax2.YTickLabel = {'-3 dB', ''};
ax2.FontSize = 12;

set(gcf, 'Position', [50 50 1300 500]);
axs = findall(gcf,'Type','axes');
set(axs, ...
    'FontName',      'Times New Roman', ...
    'FontSize',      18, ...
    'LineWidth',     1, ...
    'TickDir',       'in', ...
    'TickLabelInterpreter', 'latex', ...
    'XMinorTick',    'on', ...
    'YMinorTick',    'on', ...
    'XScale',        'log', ...
    'Layer',         'top');


set(ax, 'XTick', [2e-3 2e-2 2e-1 2]);
set(ax, 'XTickLabel', {'$2 \times 10^{-3}$','$2\times 10^{-2}$','$0.2$','2'});
