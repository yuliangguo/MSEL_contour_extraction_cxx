clear all;
addpath '/home/guoy/Desktop/Contour_Evaluation_Public/util'
% src_path = 'ALL_Extracting_Systems_Example/';
% src_path = '/home/guoy/Desktop/BSR_CFGD/Pb_TO_KGS_BSDS300/';
% folders = dir([src_path '*_1.cem']);
% src_path = '/home/guoy/Desktop/BSR_CFGD/Pb_KGS_BSDS300/';
% folders = dir([src_path '*.cem']);
% src_path = '/home/guoy/Desktop/BSR_CFGD/Pb_TO_Kovesi_BSDS300/';
% folders = dir([src_path '*.cem']);
% src_path = '/home/guoy/Desktop/BSR_CFGD/Pb_TO_SEL_BSDS300/';
% folders = dir([src_path '*.cem']);
% src_path = '/home/guoy/Desktop/BSR_CFGD/TO_SEL_BSDS300/';
% src_path = '/home/guoy/Desktop/BSR_CFGD/TO_KGS_BSDS300/';
src_path = '/home/guoy/Desktop/BSR_CFGD/TO_Kovesi_BSDS300/';
folders = dir([src_path '*.cem']);
varying_prune_len_th = [0, 5, 10, 15, 20, 30, 40, 50 ,60, 80, 100, 200, 300, 400, 500]; % varying threshold as the prune lenght of curve fragment
for l = 1:length(varying_prune_len_th)
    len_th = varying_prune_len_th(l);
    for i = 1:length(folders)
        name=[folders(i).name];
        src=[src_path,name];
        d=load_contours(src);
        I = zeros(d{1}(2),d{1}(1));
        imshow(I,'border', 'tight');
        draw_bry_contours_prune_length(d{2}, len_th, 0, [0.99,0.99,0.99]);
        h = figure(1);
    %     print_pdf([src(1:end-3), 'pdf'], h);
        set(h, 'PaperPositionMode','auto')
        print(h,'-dpng','-r0',[src(1:end-4), '_len_',num2str(len_th),'.png']);
    end

end