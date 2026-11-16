Di=1;
dat=readmatrix('Re1000_control_R75e-7_data');
Xdat2 = dat(:,2);
Ydat2 = dat(:,3);
% Load training data
tick=1;
dat(Ydat2 >= Di*1.5,:) =[];
Xdat2(Ydat2 >= Di*1.5,:) =[];
Ydat2(Ydat2 >= Di*1.5,:) =[];
dat(Ydat2 <= -Di*1.5,:) =[];
Xdat2(Ydat2 <= -Di*1.5,:) =[];
Ydat2(Ydat2 <= -Di*1.5,:) =[];
dat(Xdat2 <= -Di*2,:) =[];
Ydat2(Xdat2 <= -Di*2,:) =[];
Xdat2(Xdat2 <= -Di*2,:) =[];
dat(Xdat2 >= Di*6,:) =[];
Ydat2(Xdat2 >= Di*6,:) =[];
Xdat2(Xdat2 >= Di*6,:) =[];
Data_control=dat;

dat=readmatrix('Re1000_base_data');
Xdat3 = dat(:,2);
Ydat3 = dat(:,3);
% Load training data
tick=1;
dat(Ydat3 >= Di*1.5,:) =[];
Xdat3(Ydat3 >= Di*1.5,:) =[];
Ydat3(Ydat3 >= Di*1.5,:) =[];
dat(Ydat3 <= -Di*1.5,:) =[];
Xdat3(Ydat3 <= -Di*1.5,:) =[];
Ydat3(Ydat3 <= -Di*1.5,:) =[];
dat(Xdat3 <= -Di*2,:) =[];
Ydat3(Xdat3 <= -Di*2,:) =[];
Xdat3(Xdat3 <= -Di*2,:) =[];
dat(Xdat3 >= Di*6,:) =[];
Ydat3(Xdat3 >= Di*6,:) =[];
Xdat3(Xdat3 >= Di*6,:) =[];
Data_base=dat;

fields = {
    Data_base(:,7)/U_inf
    Data_base(:,12)/U_inf
    Data_control(:,7)/U_inf
    Data_control(:,12)/U_inf
};

fig = plotContourQuad( ...
    fields, Xdat3, Ydat3, Xdat2, Ydat2, ...
    'XLim', [-2 6], ...
    'YLim', [-1.5 1.5], ...
    'CLim', [0 2], ...
    'GridSize', [1250 400], ...
    'CylinderCentre', [0 0], ...
    'CylinderRadius', 0.5, ...
    'RecircLengths', [NaN 1.03-0.5 NaN 2.12-0.5], ...
    'RecircLabels', {'', '$~~~~~~~~~~l_c$= 1.03', '', '$l_c$= 2.12'}, ...
    'ColourBarLabel', '$u_x/U_\infty$', ...
    'PanelTitles', {'Baseline snapshot','Baseline average','Control snapshot','Control average'}, ...
    'ExportFile', 'contour_quad.pdf');

function fig = plotContourQuad(DAT, Xdat, Ydat, Xdat2, Ydat2, varargin)

p = inputParser;

p.addRequired('DAT');
p.addRequired('Xdat');
p.addRequired('Ydat');
p.addRequired('Xdat2');
p.addRequired('Ydat2');

p.addParameter('XLim', [-2 6], @(x) isnumeric(x) && numel(x) == 2);
p.addParameter('YLim', [-1.5 1.5], @(x) isnumeric(x) && numel(x) == 2);
p.addParameter('GridSize', [1200 450], @(x) isnumeric(x) && numel(x) == 2 && all(x > 1));
p.addParameter('CLim', [0 2], @(x) isnumeric(x) && numel(x) == 2);
p.addParameter('NumLevels', 40, @(x) isnumeric(x) && isscalar(x) && x >= 2);
p.addParameter('ContourLevels', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('CylinderCentre', [0 0], @(x) isnumeric(x) && numel(x) == 2);
p.addParameter('CylinderRadius', 0.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('ColourMap', redblue, @(x) isnumeric(x) && size(x,2) == 3);
p.addParameter('ColourBarLabel', '$|u|/u_\infty$', @(x) ischar(x) || isstring(x));
p.addParameter('PanelLabels', {'(a)','(b)','(c)','(d)'}, @(x) iscell(x) && numel(x) == 4);
p.addParameter('PanelTitles', {'','','',''}, @(x) iscell(x) && numel(x) == 4);
p.addParameter('FigureSize', [50 50 1150 600], @(x) isnumeric(x) && numel(x) == 4);
p.addParameter('ExportFile', '', @(x) ischar(x) || isstring(x));


p.addParameter('RecircLengths', [NaN NaN NaN NaN], @(x) isnumeric(x) && numel(x) == 4);
p.addParameter('RecircLabels', {'','','',''}, @(x) iscell(x) && numel(x) == 4);
p.addParameter('RecircLengthScale', 1.0, @(x) isnumeric(x) && isscalar(x) && x > 0);

p.parse(DAT, Xdat, Ydat, Xdat2, Ydat2, varargin{:});
opt = p.Results;


xDataCell = {Xdat(:), Xdat(:), Xdat2(:), Xdat2(:)};
yDataCell = {Ydat(:), Ydat(:), Ydat2(:), Ydat2(:)};


nx = opt.GridSize(1);
ny = opt.GridSize(2);

xGrid = linspace(opt.XLim(1), opt.XLim(2), nx);
yGrid = linspace(opt.YLim(1), opt.YLim(2), ny);
[XGrid, YGrid] = meshgrid(xGrid, yGrid);

interpolatedFields = cell(1,4);
for k = 1:4
    interpolatedFields{k} = griddata( ...
        xDataCell{k}, yDataCell{k}, DAT{k}, ...
        XGrid, YGrid, 'natural');
end

fig = figure( ...
    'Color', 'w', ...
    'Units', 'pixels', ...
    'Position', opt.FigureSize);

t = tiledlayout(fig, 2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

axesHandles = gobjects(1,4);

for k = 1:4
    ax = nexttile(t, k);
    axesHandles(k) = ax;

    hold(ax, 'on');

    contourf(ax, XGrid, YGrid, interpolatedFields{k}, opt.NumLevels, ...
        'LineColor', 'none');

    if isempty(opt.ContourLevels)
        contourLevels = linspace(opt.CLim(1), opt.CLim(2), 6);
    else
        contourLevels = opt.ContourLevels;
    end

    contour(ax, XGrid, YGrid, interpolatedFields{k}, contourLevels, ...
        'LineColor', [0.20 0.20 0.20], ...
        'LineWidth', 0.35, 'LineStyle', '-.');

    plotCylinder(ax, opt.CylinderCentre, opt.CylinderRadius);

    if isfinite(opt.RecircLengths(k))
        plotRecircLength(ax, opt.CylinderCentre, opt.CylinderRadius, ...
            opt.RecircLengths(k) * opt.RecircLengthScale, ...
            opt.RecircLabels{k});
    end

    text(ax, 0, 1.15, opt.PanelLabels{k}, ...
        'Units', 'normalized', ...
        'Interpreter', 'latex', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 14, ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'top', ...
        'Color', 'k');

    if ~isempty(opt.PanelTitles{k})
        title(ax, opt.PanelTitles{k}, ...
            'Interpreter', 'latex', ...
            'FontWeight', 'normal', ...
            'FontSize', 18);
    end

    xlim(ax, opt.XLim);
    ylim(ax, opt.YLim);
    clim(ax, opt.CLim);
    axis(ax, 'equal');

    set(ax, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 16, ...
        'LineWidth', 1.0, ...
        'TickDir', 'in', ...
        'TickLength', [0.012 0.012], ...
        'TickLabelInterpreter', 'latex', ...
        'Layer', 'top', ...
        'Box', 'on', ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on');
end

colormap(fig, opt.ColourMap);

set(axesHandles(1), 'XTickLabel', []);
set(axesHandles(2), 'XTickLabel', [], 'YTickLabel', []);
set(axesHandles(4), 'YTickLabel', []);

xlabel(t, '$x/D$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel(t, '$y/D$', 'Interpreter', 'latex', 'FontSize', 16);

cb = colorbar(axesHandles(1), 'northoutside');
cb.Layout.Tile = 'north';
cb.Layout.TileSpan = [1 2];
cb.Label.String = opt.ColourBarLabel;
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 18;
cb.Label.Rotation = 0;
cb.FontName = 'Times New Roman';
cb.FontSize = 18;
cb.TickLabelInterpreter = 'latex';

end

function plotCylinder(ax, centre, radius)
theta = linspace(0, 2*pi, 300);
xCylinder = centre(1) + radius*cos(theta);
yCylinder = centre(2) + radius*sin(theta);

fill(ax, xCylinder, yCylinder, [0.32 0.32 0.32], ...
    'EdgeColor', 'k', ...
    'LineWidth', 1.2, ...
    'FaceAlpha', 1);
end

function plotRecircLength(ax, centre, radius, lengthVal, labelStr)

x0 = centre(1) + radius;
x1 = x0 + lengthVal;
y0 = centre(2);

plot(ax, [x0 x1], [y0 y0], 'k-', 'LineWidth', 1.2);
plot(ax, [x1 x1], [y0-0.05 y0+0.05], 'k-', 'LineWidth', 0.8);

if nargin >= 5 && ~isempty(labelStr)
    text(ax, mean([x0 x1]), y0, labelStr, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 14, ...
        'FontWeight', 'normal', ...
        'Color', 'r', ...
        'Interpreter', 'latex');
end
end

exportgraphics(gcf,"figures\Fig14.pdf","ContentType","vector")