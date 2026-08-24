%load data
close all;
Re500_drag_base=readmatrix([datpath,'Re500_drag_base.txt']);
Re500_lift_base=readmatrix([datpath,'Re500_lift_base.txt']);
Re500_drag_control=readmatrix([datpath,'Re500_drag_control.txt']);
Re500_lift_control=readmatrix([datpath,'Re500_lift_control.txt']);
Re750_drag_base=readmatrix([datpath,'Re750_drag_base.txt']);
Re750_lift_base=readmatrix([datpath,'Re750_lift_base.txt']);
Re750_drag_control=readmatrix([datpath,'Re750_drag_control.txt']);
Re750_lift_control=readmatrix([datpath,'Re750_lift_control.txt']);
Re1250_drag_base=readmatrix([datpath,'Re1250_drag_base.txt']);
Re1250_lift_base=readmatrix([datpath,'Re1250_lift_base.txt']);
Re1250_drag_control=readmatrix([datpath,'Re1250_drag_control.txt']);
Re1250_lift_control=readmatrix([datpath,'Re1250_lift_control.txt']);
Re2000_drag_base=readmatrix([datpath,'Re2000_drag_base.txt']);
Re2000_lift_base=readmatrix([datpath,'Re2000_lift_base.txt']);
Re2000_drag_control=readmatrix([datpath,'Re2000_drag_control.txt']);
Re2000_lift_control=readmatrix([datpath,'Re2000_lift_control.txt']);
Re500_min_drag=(2*1.3554155e-05)/(0.5*1.225*(U_inf/2)^2); %from half sim
Re750_min_drag=(2*2.928796e-05)/(0.5*1.225*(U_inf*0.75)^2);
Re1250_min_drag=(2*7.8199976e-05)/(0.5*1.225*(U_inf*1.25)^2);
Re2000_min_drag=(2*0.00019501642)/(0.5*1.225*(U_inf*2)^2);

%process
Re500_drag_base = mean(Re500_drag_base(end-4000:end,2))/(0.5*1.225*(U_inf/2)^2);
Re750_drag_base = mean(Re750_drag_base(end-4000:end,2))/(0.5*1.225*(U_inf*0.75)^2);
Re1250_drag_base = mean(Re1250_drag_base(end-4000:end,2))/(0.5*1.225*(U_inf*1.25)^2);
Re2000_drag_base = mean(Re2000_drag_base(end-4000:end,2))/(0.5*1.225*(U_inf*2)^2);
Re500_drag_control = mean(Re500_drag_control(end-4000:end,2))/(0.5*1.225*(U_inf/2)^2);
Re750_drag_control = mean(Re750_drag_control(end-4000:end,2))/(0.5*1.225*(U_inf*0.75)^2);
Re1250_drag_control = mean(Re1250_drag_control(end-4000:end,2))/(0.5*1.225*(U_inf*1.25)^2);
Re2000_drag_control = mean(Re2000_drag_control(end-4000:end,2))/(0.5*1.225*(U_inf*2)^2);

Eta_500=Eta_calc(Re500_drag_base,Re500_drag_control,Re500_min_drag);
Eta_750=Eta_calc(Re750_drag_base,Re750_drag_control,Re750_min_drag);
Eta_1250=Eta_calc(Re1250_drag_base,Re1250_drag_control,Re1250_min_drag);
Eta_2000=Eta_calc(Re2000_drag_base,Re2000_drag_control,Re2000_min_drag);

PSD_Re{1} = (Re500_lift_base(end-4000:end,2))/(0.5*1.225*(U_inf/2)^2);
PSD_Re{2} = (Re750_lift_base(end-4000:end,2))/(0.5*1.225*(U_inf*0.75)^2);
PSD_Re{3} = (Base_lift(:,2))/(0.5*1.225*(U_inf)^2);
PSD_Re{4} = (Re1250_lift_base(end-4000:end,2))/(0.5*1.225*(U_inf*1.25)^2);
PSD_Re{5} = (Re2000_lift_base(end-4000:end,2))/(0.5*1.225*(U_inf*2)^2);
PSD_Re{6} = (Re500_lift_control(end-4000:end,2))/(0.5*1.225*(U_inf/2)^2);
PSD_Re{7} = (Re750_lift_control(end-4000:end,2))/(0.5*1.225*(U_inf*0.75)^2);
PSD_Re{8} = (R75e_7(:,end-4))/(0.5*1.225*(U_inf)^2);
PSD_Re{9} = (Re1250_lift_control(end-4000:end,2))/(0.5*1.225*(U_inf*1.25)^2);
PSD_Re{10} = (Re2000_lift_control(end-4000:end,2))/(0.5*1.225*(U_inf*2)^2);

fs = 1/0.4;
U_inf_Re = [U_inf/2 U_inf*0.75 U_inf U_inf*1.25 U_inf*2 U_inf/2 U_inf*0.75 U_inf U_inf*1.25 U_inf*2];
for i = 1:10
    n = length(PSD_Re{i});
    winLen  = max(1024, floor(0.75*n));
    window  = hann(winLen);
    overlap = round(0.5*winLen);
    NFFT    = 2^nextpow2(max(winLen, n));
    [pxx{i}, f{i}] = pwelch(PSD_Re{i}, window, overlap, NFFT, fs);
end
Reynolds = [500 750 1000 1250 2000];

Base = [Re500_drag_base Re750_drag_base 1.193797556738390 Re1250_drag_base Re2000_drag_base];
Control = [Re500_drag_control Re750_drag_control 0.883548697748480 Re1250_drag_control Re2000_drag_control];
Min = [Re500_min_drag Re750_min_drag 0.778249055510211 Re1250_min_drag Re2000_min_drag];

%plot
fig = figure('Color','w', ...
             'Units','pixels', ...
             'Position',[50 50 1250 500], ...
             'Renderer','painters');

t = tiledlayout(fig, 2, 4, ...
    'TileSpacing','compact', ...
    'Padding','compact');
ax1 = nexttile(t,[2,2]);
plot(Reynolds,Base,'-ko',LineWidth=1.5,MarkerSize=10); hold on;
plot(Reynolds,Control,'-ks',LineWidth=1.5,MarkerSize=10);
plot(Reynolds,Min,'-kx',LineWidth=1.5,MarkerSize=10); 
plot([1000 1000],[0.7 1.2],'-.r')
xlabel('Reynolds number, $Re$','Interpreter','latex')
ylabel('$\bar{C}_\mathrm{d}$','Interpreter','latex')
legend('$\bar{C}_\mathrm{d,base}$','$\bar{C}_\mathrm{d,c}$','$\bar{C}_\mathrm{d,min}$','Location','best',Interpreter='latex',fontsize=14)
grid on

tilePos = [3 4 7 8];

ax2 = nexttile(t,3);
hold(ax2,'on')
plot(f{1}/f_shed, 20*log10(pxx{1}),'-k','LineWidth', 1.5);
plot(f{6}/f_shed, 20*log10(pxx{6}),'-.r', 'LineWidth', 1.5);
grid(ax2,'on')
box(ax2,'off')
ylabel('PSD/F$^+$',Interpreter='latex')
xlim([0 15])
ylim([-175 60])
title('$Re=500$','Interpreter','latex')
legend('Base','Feedback',fontsize=12)
ax3 = nexttile(t,4);
hold(ax3,'on')
plot(f{2}/f_shed, 20*log10(pxx{2}),'-k','LineWidth', 1.5);
plot(f{7}/f_shed, 20*log10(pxx{7}),'-.r', 'LineWidth', 1.5);
grid(ax3,'on')
box(ax3,'off')
xlim([0 15])
ylim([-175 60])
title('$Re=750$','Interpreter','latex')
ax4 = nexttile(t,7);
hold(ax4,'on')
plot(f{4}/f_shed, 20*log10(pxx{4}),'-k','LineWidth', 1.5);
plot(f{9}/f_shed, 20*log10(pxx{9}),'-.r', 'LineWidth', 1.5);
grid(ax4,'on')
box(ax4,'off')
xlabel('F$^+$',Interpreter='latex')
ylabel('PSD/F$^+$',Interpreter='latex')
xlim([0 15])
ylim([-175 60])
title('$Re=1250$','Interpreter','latex')
ax5 = nexttile(t,8);
hold(ax5,'on')
plot(f{5}/f_shed, 20*log10(pxx{5}),'-k','LineWidth', 1.5);
plot(f{10}/f_shed, 20*log10(pxx{10}),'-.r', 'LineWidth', 1.5);
grid(ax5,'on')
box(ax5,'off')
xlabel('F$^+$',Interpreter='latex')
xlim([0 15])
ylim([-175 60])
title('$Re=2000$','Interpreter','latex')
axs = [ax1 ax2 ax3 ax4 ax5];
set(axs, ...
    'FontName',      'Times New Roman', ...
    'FontSize',      18, ...
    'LineWidth',     1, ...
    'TickDir',       'in', ...
    'TickLabelInterpreter', 'latex', ...
    'XMinorTick',    'on', ...
    'YMinorTick',    'on', ...
    'XScale',        'linear', ...
    'Layer',         'bottom', ...
    'Box','off');

exportgraphics(gcf,"Fig17.eps","ContentType","vector")