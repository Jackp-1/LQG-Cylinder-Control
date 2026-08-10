fig = figure('Color','w', ...
             'Units','pixels', ...
             'Position',[50 50 1250 400], ...
             'Renderer','painters');

t = tiledlayout(fig,1,4, ...
    'TileSpacing','compact', ...
    'Padding','compact');
t.Units = 'normalized';
t.OuterPosition = [0.03 0.02 0.97 0.96];
lw  = 1.6;
clr1 = [0 0 0];
clr2 = [230/255 20/255   0/255];
clr3 = [0/255   200/255   0/255];

ax1 = nexttile(t,1);
yyaxis left
ax = gca;
ax.YColor = 'k';
semilogy(1:length(Sd),Sd,'-', 'Color', clr1, 'LineWidth', lw);
hold on;
scatter(1:length(Sd),Sd,30,'black');
scatter(1:r,Sd(1:r),30,'red')
set(gca,'yscale','log')
ylabel('Singular value, $\sigma_r$','Interpreter','latex',Rotation=90)
ylim([0.1 30])
grid(ax, "on")
yyaxis right
ax = gca;
ax.YColor = [0/255 20/255   230/255];
semilogy(cumsum(Sd)/sum(Sd),'--','Color',[0/255 0/255 0/255], 'LineWidth', lw)
hold on;
scatter(1:length(Sd),cumsum(Sd/sum(Sd)),30,'s','MarkerEdgeColor',[0/255 0/255 0/255]);
scatter(1:r,cumsum(Sd(1:r)/sum(Sd)),30,'s','MarkerEdgeColor',[230/255 20/255 0/255])
plot(11*ones(1,10),linspace(0,0.7723,10),'-.',LineWidth=1.5); plot(linspace(11,25,10),ones(1,10)*0.7723,'-.',LineWidth=1.5);
xlim([0 25])
ylabel('Cumulative Energy','Interpreter','latex',Rotation=270)
grid(ax, 'on');
box(ax1, 'off');
xlabel('r','Interpreter','latex')
ylim([0 1])

ax2 = nexttile(t,2,[1 3]);
nx=2400; %Number of grid points in x direction
ny=900; %Number of grid points in y direction
D=1;
X = linspace(-2*D,6*D,nx);
Y = linspace(-1.5*D,1.5*D,ny);
[X,Y]=meshgrid(X,Y);
[X,Y,DAT]=griddata(Xdat,Ydat,Phi{r}(:,7),X,Y);
DAT(DAT == 0) = NaN;
r=11;
contourf(real(DAT),40,'LineColor','none'); hold on; 
contour(real(DAT),6,'LineColor','k','LineStyle','-.')
colormap('redblue')
title('Mode 1 (25.3\%), St=0.21','FontWeight','Normal','FontSize',16)
hold on

% Clean up axes
set(gca,'XTick',[25 300 600 900 1200 1500 1800 2100 2375],'XTickLabel',{'-2','-1','0','1','2','3','4','5','6'});
set(gca,'YTick',[150 450 750],'YTickLabel',{'1','0','-1'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times')
xlabel('$x/ D$','Interpreter','latex')
ylabel('$y/ D$','rotation',90,'Interpreter','latex')
%Plot Cylinder overlay
theta= (1:101)/100'*2*pi;
x=600+150*sin(theta);
y=450+150*cos(theta);
fill(x,y,[.3 .3 .3],'FaceAlpha',1)% place cylinder
plot(x,y,'k','LineWidth',1.2)% cylinder boundary
xlim([25 2375])
ylim([5 895])
axs = findall(gcf,'Type','axes');
set(axs, ...
    'FontName',      'Times New Roman', ...
    'FontSize',      18, ...
    'LineWidth',     1, ...
    'TickDir',       'in', ...
    'TickLabelInterpreter', 'latex', ...
    'XMinorTick',    'on', ...
    'YMinorTick',    'on', ...
    'XScale',        'linear', ...
    'Layer',         'bottom');

set(ax2,'Layer','top')

exportgraphics(gcf,"Fig4.eps","ContentType","vector")