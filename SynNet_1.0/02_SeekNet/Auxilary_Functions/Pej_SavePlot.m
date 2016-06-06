% This scripts saves the figure, Fig, in the address, OutFig, in .eps, and
% .fig formats and closes it.

function Pej_SavePlot(Fig, OutFig, Transparent)
if nargin <3; Transparent = false; end

[pathstr,~,~] = fileparts(OutFig);
if ~isempty(pathstr) && ~exist(pathstr, 'dir'); mkdir(pathstr);end

saveas(Fig, [OutFig '.fig'])
if Transparent
    set(Fig, 'color', 'none',         'inverthardcopy', 'off');
end
% try
% set(Fig, 'PaperPositionMode', 'auto');
% catch err
% end
% print(Fig, '-depsc2', '-r0', [OutFig '.eps']);


set(Fig, 'PaperUnits','centimeters');
set(Fig, 'Units','centimeters');
pos=get(Fig,'Position');
set(Fig, 'PaperSize', [pos(3) pos(4)]*1.05);
set(Fig, 'PaperPositionMode', 'manual');
set(Fig, 'PaperPosition',[.025*pos(3) .025*pos(4) pos(3) pos(4)]);

print(Fig, '-dpdf', [OutFig '.pdf']);
close(Fig)
disp(['Figure saved: ' OutFig] )
end