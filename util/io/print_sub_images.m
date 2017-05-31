clear all;
src_path = '/home/guoy/Desktop/contour_shi_new/divide_merge_test/';
files = dir([src_path '*png']);
n = 60;
for i = 1:length(files)
    i
    name=[files(i).name];
    cem_src=[src_path,name(1:end-3),'cem'];
    d=load_contours(cem_src);
%     [edg, edgemap, thetamap] = load_edg([src(1:end-4) '_1.edg']);
%     I = zeros(d{1}(2),d{1}(1),3);
    I = imread([src_path, name]);
    imshow(I,'border', 'tight');

    draw_contours_first_n(d{2}, 100);
%     disp_edg(edg);
    h = figure(1);
%     keyboard;
%     print_pdf([src(1:end-3), 'pdf'], h);
    set(h, 'PaperPositionMode','auto')
    print(h,'-djpeg95','-r0',[cem_src(1:end-4), '_sub_c.jpg']);
end