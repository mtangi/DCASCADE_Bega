function map = shadesOfColor(color, ncolors)
% shadesOfColor of color generates a fading colorscale from variable
% "color"

HSV = rgb2hsv(color);
shades = log(linspace(exp(.2),exp(0.9),ncolors));
%shades = linspace(.11,.88,ncolors);
colors = zeros(ncolors, 3);
colors(:,1) = HSV(1);
colors(:,2) = shades;
colors(:,3) = HSV(3);
map = hsv2rgb(colors);
end