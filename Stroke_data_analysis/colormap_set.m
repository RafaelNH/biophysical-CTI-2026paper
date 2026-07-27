colormap(mycmap)
pos = get(ax, 'position');
cb = colorbar('southoutside');
cb.Position =  [pos(1) pos(2) pos(3)+0.005 0.01];
cb.FontSize = 8;