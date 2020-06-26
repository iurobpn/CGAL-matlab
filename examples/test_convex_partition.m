%% Example of usage of polygon convex partition functions of CGAL library
% in Matlab (Tested with Matlab R2019b)
%
%   Iuro Nascimento 25/06/2020
%
close all
clear all
poly = [25,150;
          40,185;
          51,173;
          59,174;
          60,176;
          50,193;
          88,198;
          81,175;
          86,169;
          92,167;
          110,185;
          125,150];

figure
subplot(3,1,1)
plot_polygon(poly);
title('original polygon')
axis equal

N = size(poly,1);
ps = poly(N:-1:1,:); %the vertices have to be in sorted in counter-clockwise order

Yopt = optimal_convex_partition(ps);
subplot(3,1,2);
plot_polygon(Yopt)
title('optimal convex partition')
axis equal

Ygreene = greene_approx_convex_partition(ps);
subplot(3,1,3);
plot_polygon(Ygreene)
title('Greene partition')
axis equal