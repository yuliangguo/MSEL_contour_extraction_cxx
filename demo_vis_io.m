addpath (genpath('util/'));

img = imread('10081.jpg');
[h,w,~]= size(img);

% load in edges and contours
[edges, edgemap, thetamap] = load_edg('10081.edg');
[CEM, edges, cfrags_idx] = load_contours('10081.cem');

% visualize edges
figure(1);
imshow(edgemap, 'border', 'tight');

% visualize contours
figure(2);
imshow(img, 'border', 'tight'); hold on;
draw_contours(CEM{2});

% output edges and contour files
write_cem('10081.cem', CEM{2}, h, w)
save_edg('10081.edg', edges, [w,h])