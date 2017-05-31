clear all; close all;
src_path = '/home/guoy/Desktop/Skeleton_Track_Project/mouse_data/Mouse_Ground_Truth/crim_contour_w/';
files = dir([src_path '*cem']);

for i = 1:length(files)
    cem_src = [src_path, files(i).name];
    d=load_contours(cem_src);
    
    I = zeros(d{1,1}(1,2), d{1,1}(1,1));
    h = figure(1);
    imshow(I, 'border', 'tight');
    hold on;
    draw_contours(d{2}, 0, 0, [0.99,0.99,0.99]); %draw white
    
    set(h, 'PaperPositionMode','auto')
    print(h,'-djpeg95','-r0',[cem_src(1:end-3), 'jpg']);
end
