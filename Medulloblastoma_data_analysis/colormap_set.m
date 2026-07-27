% colormap(mycmap)
% pos = get(ax, 'position');
% cb = colorbar('southoutside');
% cb.Position =  [pos(1) pos(2)+0.005 pos(3) 0.01];
% cb.FontSize = fff;

colormap(mycmap)
pos = get(ax, 'position');
cb = colorbar('eastoutside');
cb.Position =  [pos(1)+pos(3)*0.92,...
    pos(2)+pos(4)*0.05,...
    0.01,...
    pos(4)*0.9];
cb.FontSize = fff*0.8;