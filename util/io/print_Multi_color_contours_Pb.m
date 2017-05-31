clear all; close all;
% pb_path = '/home/guoy/Desktop/contour_shi_new/BSDS300_test/Pb_base/';
% img_path = '/home/guoy/Desktop/BSR_CFGD/BSDS300_test_img/';
% src_path = '/home/guoy/Desktop/BSR_CFGD/Pb_KGS_BSDS300/';
pb_path = '/home/guoy/Desktop/Print_PDF/Pb_KGS_CFGD/';
img_path = '/home/guoy/Desktop/contour_shi_new/CFGD_test_1/';
src_path = '/home/guoy/Desktop/Print_PDF/Pb_KGS_CFGD/';
folders = dir([src_path, '*cem']);
pb_th = 0.1250;
for i = 1:length(folders)
    src= [src_path folders(i).name];
    d=load_contours(src);
%     I = zeros(d{1}(2),d{1}(1),3);
    pb_file_name = [pb_path, folders(i).name(1:end-3) 'mat'];
    load(pb_file_name);
    pb = pb/max(max(pb));
    
    I=imread([img_path folders(i).name(1:end-3) 'jpg']);
    imshow(I,'border', 'tight');
%     draw_contours(d{2});
    draw_multicolor_contours_prune_pb(d{2,1}, pb, pb_th);
    h = figure(1);
    set(h, 'PaperPositionMode','auto')
    print_pdf([src(1:end-3), 'pdf']);
%     print(h,'-djpeg95','-r0',[src(1:end-4), '_color.jpg']);
end