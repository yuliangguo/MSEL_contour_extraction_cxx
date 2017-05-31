clear all; close all;
% src_path = '/media/NewVolume_1/contour_shi_new/CFGD_test_1/PB_TO_base/Refined_Results/';
% % src_path = '/media/NewVolume_1/contour_shi_new/CFGD_test_1/Pb_TO_base/';
% edge_path = '/media/NewVolume_1/contour_shi_new/CFGD_test_1/PB_TO_base/';
% image_path = '/media/NewVolume_1/contour_shi_new/CFGD_test_1/';
% 
% files = dir([image_path '*jpg']);
% for i = 1:length(files)
%     i
%     name=[files(i).name];
%     src=[src_path,name(1:end-4), '_refined.cem'];
%     d=load_contours(src);
%     [edg, edgemap, thetamap] = load_edg([edge_path, name(1:end-4), '_1.edg']);
% %     I = zeros(d{1}(2),d{1}(1),3);
%     I = imread([image_path, name]);
%     imshow(I,'border', 'tight');
% 
%     draw_contours(d{2});
%     disp_edg(edg);
%     h = figure(1);
% %     keyboard;
%     print_pdf([src(1:end-4), '_exact.pdf'], h);
% %     set(h, 'PaperPositionMode','auto')
% %     print(h,'-djpeg95','-r0',[src(1:end-4), '_color.jpg']);
% end

clear all; close all;
addpath ../Print_PDF/
src_path = '../CFGD/Coarse_Scale/curve_fragments_oriented/';
% src_path = '/media/NewVolume_1/contour_shi_new/CFGD_test_1/Pb_TO_base/';
% edge_path = 'CFGD_test_1/TO_base/frags_hyp_edgs/';
image_path = '../CFGD/Coarse_Scale/GT_img/';

files = dir([image_path '*jpg']);
for i = 1:length(files)
    i
    name=[files(i).name];
    src=[src_path,name(1:end-4), '_s1.cem'];
    d=load_contours(src);
%     load([src(1:end-4) '.mat']);
%     [edg, edgemap, thetamap] = load_edg([edge_path, name(1:end-4), '.edg']);
%     I = zeros(d{1}(2),d{1}(1),3);
    I = imread([image_path, name]);
    imshow(I,'border', 'tight');

    draw_contours(d{2});
%     disp_edg(edg);
    h = figure(1);
%     keyboard;
    print_pdf([src(1:end-4), '_s1.pdf'], h);
%     set(h, 'PaperPositionMode','auto')
%     print(h,'-djpeg95','-r0',[src(1:end-4), '_color.jpg']);
end