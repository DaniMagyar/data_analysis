function [colors] = BAfc_colors

colors.PN_primary = [8,104,172]/255;
colors.PN_secondary = [43,140,190]/255;
colors.IN_primary = [174,1,126]/255;
colors.IN_secondary = [247,104,161]/255;

colors.raw_primary = [44,162,95]/255;
colors.raw_secondary = [153,216,201]/255;


% Colormap blue:white:red -------------------------------------------------
c1 = 1/255*[0,115,185];
c2 = 1/255*[239,239,239];
c3 = 1/255*[217,84,26];
mycolormap(1:20,1)  = linspace(c1(1),c2(1),20);
mycolormap(21:64,1) = linspace(c2(1),c3(1),44);
mycolormap(1:20,2)  = linspace(c1(2),c2(2),20);
mycolormap(21:64,2) = linspace(c2(2),c3(2),44);
mycolormap(1:20,3)  = linspace(c1(3),c2(3),20);
mycolormap(21:64,3) = linspace(c2(3),c3(3),44);  
colors.Heatmap = mycolormap;
colors.c1 = c1;
colors.c3 = c3;