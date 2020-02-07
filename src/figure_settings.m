function [FONTSIZE, LINEWIDTH, FIGSIZE, SCATTER ] = plot_settings( fsize, lw, figsize, scatter )
%plot_settings : plot_settings(FONTSIZE, LINEWIDTH, FIGSIZE, SCATTER)
%   Detailed explanation goes here

set(groot, 'DefaultTextFontSize', fsize);
set(groot, 'DefaultAxesFontSize', fsize);
% set(0, 'DefaultAxesFontName', 'Verdana');
set(groot, 'DefaultAxesLineWidth', lw);
set(groot, 'DefaultAxesTickLength', [0.02 0.025]);
set(groot, 'DefaultAxesBox', 'on');
set(groot, 'DefaultLineLineWidth', lw);
set(groot, 'DefaultLineMarkerSize', 5);
set(groot, 'DefaultPatchLineWidth', lw);
set(groot, 'DefaultFigureColor', [1 1 1]);
set(groot, 'DefaultFigurePaperPosition',[2.65 4.3 3.2 2.4]);
set(groot, 'DefaultFigureDockControls', 'off')
set(groot, 'DefaultLineMarkerSize', 2);
set(groot, 'DefaultLegendBox', 'off');
set(groot, 'DefaultAxesColorOrder',[[0 0 0]; [0 0 0]])

% scatter properties
set(groot, 'DefaultScatterLineWidth', lw);

% SET FIGURE POSITION
set(groot, 'DefaultFigurePosition', [800 300 figsize])

FONTSIZE = fsize;
LINEWIDTH = lw;
FIGSIZE = figsize;
SCATTER = scatter;
end

