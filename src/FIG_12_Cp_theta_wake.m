x1_avg_base=sortrows(readmatrix([datpath, 'x1_u_avg_base.csv']),2);
x1_avg_base(:,1) = x1_avg_base(:,1)/U_inf;
x15_avg_base=sortrows(readmatrix([datpath, 'x15_u_avg_base.csv']),2);
x15_avg_base(:,1) = x15_avg_base(:,1)/U_inf;
x2_avg_base=sortrows(readmatrix([datpath, 'x2_u_avg_base.csv']),2);
x2_avg_base(:,1) = x2_avg_base(:,1)/U_inf;
x25_avg_base=sortrows(readmatrix([datpath, 'x25_u_avg_base.csv']),2);
x25_avg_base(:,1) = x25_avg_base(:,1)/U_inf;
x3_avg_base=sortrows(readmatrix([datpath, 'x3_u_avg_base.csv']),2);
x3_avg_base(:,1) = x3_avg_base(:,1)/U_inf;
x1_avg_control=sortrows(readmatrix([datpath, 'x1_u_avg_control.csv']),2);
x1_avg_control(:,1) = x1_avg_control(:,1)/U_inf;
x15_avg_control=sortrows(readmatrix([datpath, 'x15_u_avg_control.csv']),2);
x15_avg_control(:,1) = x15_avg_control(:,1)/U_inf;
x2_avg_control=sortrows(readmatrix([datpath, 'x2_u_avg_control.csv']),2);
x2_avg_control(:,1) = x2_avg_control(:,1)/U_inf;
x25_avg_control=sortrows(readmatrix([datpath, 'x25_u_avg_control.csv']),2);
x25_avg_control(:,1) = x25_avg_control(:,1)/U_inf;
x3_avg_control=sortrows(readmatrix([datpath, 'x3_u_avg_control.csv']),2);
x3_avg_control(:,1) = x3_avg_control(:,1)/U_inf;
Theta_Cp_Base=readmatrix([datpath, 'Theta_Cp_Base.csv']);
Theta_Cp_Base = sortrows(Theta_Cp_Base,1);
Theta_Cp_Control=readmatrix([datpath, 'Theta_Cp_Control.csv']);
Theta_Cp_Control = sortrows(Theta_Cp_Control,1);


fig = figure('Color','w', ...
             'Units','pixels', ...
             'Position',[50 50 1250 400], ...
             'Renderer','painters');

t = tiledlayout(fig, 1, 9, ...
    'TileSpacing','compact', ...
    'Padding','compact');

ax1 = nexttile(t,1,[1 4]);
axis equal
plot(Theta_Cp_Base(1:264,1)/180,Theta_Cp_Base(1:264,2),'-k','LineWidth',1.5); hold on;
plot(Theta_Cp_Base(265:496,1)/180,Theta_Cp_Base(265:496,2),'-k','LineWidth',1.5);
plot(Theta_Cp_Base(497:end,1)/180,Theta_Cp_Base(497:end,2),'-k','LineWidth',1.5);
plot(Theta_Cp_Control(1:264,1)/180,Theta_Cp_Control(1:264,2),'-.r','LineWidth',1.5)
plot(Theta_Cp_Control(265:496,1)/180,Theta_Cp_Control(265:496,2),'-.r','LineWidth',1.5); 
plot(Theta_Cp_Control(497:end,1)/180,Theta_Cp_Control(497:end,2),'-.r','LineWidth',1.5); 
patch(ax1, ...
    [-0.392565 -0.294963 -0.294963 -0.392565], ...
    [-1.5 -1.5 1.5 1.5], ...
    [1 0 0], ...
    'FaceAlpha',0.18, ...
    'EdgeColor','none', ...
    'HandleVisibility','off');
patch(ax1, ...
    [0.392565 0.294963 0.294963 0.392565], ...
    [-1.5 -1.5 1.5 1.5], ...
    [1 0 0], ...
    'FaceAlpha',0.18, ...
    'EdgeColor','none', ...
    'HandleVisibility','off');
text(ax1,-0.294963-.03,0,'Actuation Slot', ...
    'Rotation',90, ...
    'Color',[0 0 0], ...
    'FontName','Times New Roman', ...
    'FontSize',14, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle');
text(ax1,0.392565-.03,0,'Actuation Slot', ...
    'Rotation',90, ...
    'Color',[0 0 0], ...
    'FontName','Times New Roman', ...
    'FontSize',14, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle');
annotation("doublearrow",[0.245,0.245],[0.2,0.33])
text(-0.23,-0.95,'$\Delta C_{p,b}~~ 44.5\%$','FontSize',14,Interpreter='latex')
legend({'Baseline','','','Control'}, ...
    'Interpreter','latex', ...
    'FontSize', 12, ...
    'Location','northwest', ...
    'Box','on');
grid on;
box(ax1, 'off');
xlabel('$\theta/\pi$','Interpreter','latex',fontsize=16)
ylabel('$C_p$','Interpreter','latex',fontsize=16)

ax2 = nexttile(t);
plot(x1_avg_base(:,1),x1_avg_base(:,2),'-k','LineWidth',1.5); hold on;
plot(x1_avg_control(:,1),x1_avg_control(:,2),'-.r','LineWidth',1.5);
xlim([-0.5 2])
ylim([0 1.5])
title('$x/D=1$','FontSize',14,'Interpreter','latex')
ylabel('$y/D$','FontSize',14,Interpreter='latex')
xlabel('$u/U_\infty$','FontSize', 14,'Interpreter','latex')
grid on;
box(ax2, 'off');

ax3 = nexttile(t);
plot(x15_avg_base(:,1),x15_avg_base(:,2),'-k','LineWidth',1.5); hold on
plot(x15_avg_control(:,1),x15_avg_control(:,2),'-.r','LineWidth',1.5);
xlim([-0.5 2])
ylim([0 1.5])
title('$x/D=1.5$','FontSize',14,'Interpreter','latex')
xlabel('$u/U_\infty$','FontSize', 14,'Interpreter','latex')
grid on;
box(ax3, 'off');

ax4 = nexttile(t);
plot(x2_avg_base(:,1),x2_avg_base(:,2),'-k','LineWidth',1.5); hold on
plot(x2_avg_control(:,1),x2_avg_control(:,2),'-.r','LineWidth',1.5);
xlim([-0.5 2])
ylim([0 1.5])
title('$x/D=2$','FontSize',14,'Interpreter','latex')
xlabel('$u/U_\infty$','FontSize', 14,'Interpreter','latex')
grid on;
box(ax4, 'off');

ax5 = nexttile(t);
plot(x25_avg_base(:,1),x25_avg_base(:,2),'-k','LineWidth',1.5); hold on
plot(x25_avg_control(:,1),x25_avg_control(:,2),'-.r','LineWidth',1.5);
xlim([-0.5 2])
ylim([0 1.5])
title('$x/D=2.5$','FontSize',14,'Interpreter','latex')
xlabel('$u/U_\infty$','FontSize', 14,'Interpreter','latex')
grid on;
box(ax5, 'off');

ax6 = nexttile(t);
plot(x3_avg_base(:,1),x3_avg_base(:,2),'-k','LineWidth',1.5); hold on
plot(x3_avg_control(:,1),x3_avg_control(:,2),'-.r','LineWidth',1.5);
xlim([-0.5 2])
ylim([0 1.5])
title('$x/D=3$','FontSize',14,'Interpreter','latex')
xlabel('$u/U_\infty$','FontSize', 14,'Interpreter','latex')
grid on;
box(ax6, 'off');

axs = [ax1 ax2 ax3 ax4 ax5 ax6];
set(axs, ...
    'FontName',      'Times New Roman', ...
    'FontSize',      12, ...
    'LineWidth',     1, ...
    'TickDir',       'in', ...
    'TickLabelInterpreter', 'latex', ...
    'XMinorTick',    'on', ...
    'YMinorTick',    'on', ...
    'XScale',        'linear', ...
    'Layer',         'top');
