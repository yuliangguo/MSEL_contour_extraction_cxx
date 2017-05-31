clear all; close all;
addpath '/home/guoy/Desktop/Contour_Evaluation_Public/util';
% pb_path = '/home/guoy/Desktop/contour_shi_new/BSDS300_test/Pb_base/';
img_path = '/home/guoy/Desktop/contour_shi_new/CFGD_test_1/';
src_path = '/home/guoy/Desktop/Print_PDF/TO_SEL_CFGD/';
folders = dir([src_path, '*cem']);
len_th = 0;
for i = 1:length(folders)
    src= [src_path folders(i).name];
    d=load_contours(src);
%     I = zeros(d{1}(2),d{1}(1),3);
%     pb_file_name = [pb_path, folders(i).name(1:end-3) 'mat'];
%     load(pb_file_name);
%     pb = pb/max(max(pb));
    
    I=imread([img_path folders(i).name(1:end-3) 'jpg']);
    imshow(I,'border', 'tight');
%     draw_contours(d{2});
    draw_multicolor_contours_prune_length(d{2,1}, len_th);
    h = figure(1);
    set(h, 'PaperPositionMode','auto')
    print_pdf([src(1:end-3), 'pdf']);
%     print(h,'-djpeg95','-r0',[src(1:end-4), '_color.jpg']);
end