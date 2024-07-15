function saveFigure(handle, fileName, width, height)

set(handle, 'paperunits','inches');
set(handle, 'papersize',[width height]);
set(handle, 'paperposition',[0 0 width height]);
saveas(handle, fileName);

end