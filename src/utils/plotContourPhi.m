function [f1, Col, Row]=plotContourPhi(DAT,Xdat,Ydat,val)
f1 = figure;
nx=2400; %Number of grid points in x direction
ny=900; %Number of grid points in y direction
D=1;
X = linspace(-2*D,6*D,nx);
Y = linspace(-1.5*D,1.5*D,ny);
[X,Y]=meshgrid(X,Y);
[X,Y,DAT]=griddata(Xdat,Ydat,DAT(:,end),X,Y);
DAT(DAT == 0) = NaN;
contourf(DAT,50,'LineColor','none');
%clim([-1.5e-4 1.5e-4]);
%cb = colorbar('northoutside'); 
colormap('redblue')
hold on
% Clean up axes
set(gca,'XTick',[25 300 600 900 1200 1500 1800 2100 2375],'XTickLabel',{'-2','-1','0','1','2','3','4','5','6'});
set(gca,'YTick',[150 450 750],'YTickLabel',{'1','0','-1'});
set(gcf,'Position',[0 0 1200 600])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)
axis equal

%contour(DAT,[0.0025 0.0075 0.0125],':k','LineWidth',1.1) %vel
%contour(DAT,[-0.0025 -0.0075 -0.0125],'-.k','LineWidth',1.1) %vel
xlabel('$x/ \phi$','fontsize',19,'Interpreter','latex')
ylabel('$y/ \phi$','fontsize',19,'rotation',0,'Interpreter','latex')
%Plot Cylinder overlay
theta= (1:101)/100'*2*pi;
x=600+150*sin(theta);
y=450+150*cos(theta);
fill(x,y,[.3 .3 .3],'FaceAlpha',1)% place cylinder
plot(x,y,'k','LineWidth',1.2)% cylinder boundary
for i=1:size(val,1)
[~,row]=min(abs(Y-Ydat(val(i))));
[~,col]=min(abs(X'-Xdat(val(i))));
Col(i) = X(1,col(1));
Row(i) = Y(row(1),1);
% Plot Sensor Location
hold on;
%scatter(col(1),row(1),50,"filled",'m','LineWidth',0.5,'MarkerEdgeColor','k')
end


xlim([25 2375])
ylim([5 895])
end