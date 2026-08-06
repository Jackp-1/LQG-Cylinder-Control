%% Plot contour
function [f1, Col, Row]=plotContour(DAT,Xdat,Ydat,val,plotTitle)
%f1 = figure;
nx=2400; %Number of grid points in x direction
ny=900; %Number of grid points in y direction
D=1;
X = linspace(-2*D,6*D,nx);
Y = linspace(-1.5*D,1.5*D,ny);
[X,Y]=meshgrid(X,Y);
[X,Y,DAT]=griddata(Xdat,Ydat,DAT(:,end),X,Y);
DAT(DAT == 0) = NaN;
contourf(DAT, 40, 'LineColor', 'none');
%cb = colorbar('northoutside'); 
%ylabel(cb, '$\vert \mathbf{u} \vert / u_\infty$','Interpreter','latex','FontSize',19,'Rotation',0) %Normal velocity
%clim([0 1.4]);
colormap('turbo')
hold on
% Clean up axes
set(gca,'XTick',[25 300 600 900 1200 1500 1800 2100 2375],'XTickLabel',{'-2','-1','0','1','2','3','4','5','6'});
set(gca,'YTick',[150 450 750],'YTickLabel',{'1','0','-1'});
set(gcf,'Position',[0 0 1200 600])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)
axis equal

contour(DAT, [0.2 0.4 0.6 0.8 1.0 1.2], ...
    'LineColor', [0.25 0.25 0.25], ...
    'LineWidth', 0.4);
xlabel('$x/ \phi$','fontsize',19,'Interpreter','latex')
ylabel('$y/ \phi$','fontsize',19,'rotation',90,'Interpreter','latex')
%Plot Cylinder overlay
theta= (1:101)/100'*2*pi;
x=600+150*sin(theta);
y=450+150*cos(theta);
fill(x,y,[.3 .3 .3],'FaceAlpha',1)% place cylinder
plot(x,y,'k','LineWidth',1.2)% cylinder boundary

% Plot Sensor Location
if val~=0
    for i=1:size(val,1)
        [~,row]=min(abs(Y-Ydat(val(i))));
        [~,col]=min(abs(X'-Xdat(val(i))));
        Col(i) = X(1,col(1));
        Row(i) = Y(row(1),1);
        %scatter(col(1),row(1),50,"filled",'m','LineWidth',0.5,'MarkerEdgeColor','k'); hold on
    end
    Col
    Row
end

title(plotTitle, 'Interpreter','latex', ...
    'FontWeight','normal', ...
    'FontSize',16);

xlim([25 2375])
ylim([5 895])
end