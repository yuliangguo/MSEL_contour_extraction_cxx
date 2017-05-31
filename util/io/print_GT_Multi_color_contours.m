clear all;
src_path = 'TO_SEL_CFGD/';
folders = dir([src_path '*.cem']);
for i = 1:length(folders)
    name= folders(i).name ;
    src=[src_path,name];
    d=load_contours(src);
%     I = zeros(d{1}(2),d{1}(1),3);
    I = imread([src_path name(1:end-4) '.jpg']);
    h = figure(1);
    imshow(I,'border', 'tight');
    draw_contours(d{2});
    print_pdf([src(1:end-4), '.pdf'], h);
%     set(h, 'PaperPositionMode','auto')
%     print(h,'-djpeg95','-r0',[src(1:end-4), '_color.jpg']);
end