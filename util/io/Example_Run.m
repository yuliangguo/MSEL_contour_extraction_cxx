clear all; close all;
I=imread('/home/guoy/Desktop/contour_extraction_revisted/figs/16068.jpg');
d=load_contours('/home/guoy/Desktop/contour_extraction_revisted/figs/16068.cem');
% edges = load_edg('data/35049.edg');
% d2=load_contours('syn_data/contImage4.cem');
% I = zeros(d2{1}(2),d2{1}(1),3);
h = figure(1);
imshow(I, 'border', 'tight');
hold on;
%draw_red_contours(d{2});
% draw_green_contours(d2{2});
% disp_edg(edges)
% draw_contours(d{2})
draw_contours(d2{2}, 0, 0, [0.99,0.99,0.99]); %draw white 
% print(h, '-dpng', 'syn_data/contImage4_gt.png'); % never use .jpg
