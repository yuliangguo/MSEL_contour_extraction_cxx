clear all; close all;
src_path = '/media/NewVolume_1/contour_shi_new/CFGD_test_1/PB_TO_base/';
% src_path = '/media/NewVolume_1/contour_shi_new/BSDS300_test/Pb_TO_base/';

files = dir([src_path '*cem']);
for i = 1:length(files)
    i
    name=[files(i).name];
    src=[src_path,name];
    d=load_contours(src);
    [edg, edgemap, thetamap] = load_edg([src(1:end-4) '_1.edg']);
%     I = zeros(d{1}(2),d{1}(1),3);
    I = imread([src_path, '../', name(1:end-3), 'jpg']);
    imshow(I,'border', 'tight');

    draw_contours(d{2});
    disp_edg(edg);
    h = figure(1);
%     keyboard;
    print_pdf([src(1:end-4), '_exact.pdf'], h);
%     set(h, 'PaperPositionMode','auto')
%     print(h,'-djpeg95','-r0',[src(1:end-4), '_color.jpg']);
end

% src_path = '/home/guoy/Desktop/Compute_contours_scripts/PB_TO_BSDS300/';
% files = dir(src_path);
% for i = 3:length(files)
%     i
%     name=[files(i).name];
%     src=[src_path,name,'/', name, '.cem'];
%     d=load_contours(src);
%     [edg, edgemap, thetamap] = load_edg([src(1:end-4) '.edg']);
% %     I = zeros(d{1}(2),d{1}(1),3);
%     I = imread([src(1:end-4), '.jpg']);
%     imshow(I,'border', 'tight');
% 
%     draw_contours(d{2});
%     disp_edg(edg);
%     h = figure(1);
% %     keyboard;
%     print_pdf([src(1:end-3), 'pdf'], h);
% %     set(h, 'PaperPositionMode','auto')
% %     print(h,'-djpeg95','-r0',[src(1:end-4), '_color.jpg']);
% end