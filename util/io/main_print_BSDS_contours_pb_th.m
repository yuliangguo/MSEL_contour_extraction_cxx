clear all;
addpath '/home/guoy/Desktop/Contour_Evaluation_Public/util'

pb_path = '/home/guoy/Desktop/contour_shi_new/BSDS300_test/Pb_TO_base/';

% src_path = 'ALL_Extracting_Systems_Example/';
% src_path = '/home/guoy/Desktop/BSR_CFGD/Pb_TO_KGS_BSDS300/';
% src_path = '/home/guoy/Desktop/contour_shi_new/BSDS300_test/Pb_TO_base/';
% folders = dir([src_path '*_1.cem']);
% src_path = '/home/guoy/Desktop/contour_shi_new/BSDS300_test/Pb_base/';
% folders = dir([src_path '*.cem']);
% src_path = '/home/guoy/Desktop/BSR_CFGD/Pb_TO_Kovesi_BSDS300/';
% folders = dir([src_path '*.cem']);
src_path = '/home/guoy/Desktop/BSR_CFGD/Pb_TO_SEL_BSDS300/';
folders = dir([src_path '*.cem']);
% nthresh = 15;
% varying_prune_len_th = [0, 3, 5, 10, 15, 20, 25, 30, 40, 50 ,60, 70, 100, 200, 500]; % varying threshold as the prune lenght of curve fragment
% thresh = linspace(1/(nthresh+1),1-1/(nthresh+1),nthresh)';
% for l = 1:length(thresh)
%     th = thresh(l);
    for i = 1:length(folders)
        i
        name=[folders(i).name]
        src=[src_path,name];
        d=load_contours(src);
        pb_file_name = [pb_path, name(1:end-3) 'mat'];
        load(pb_file_name);
        I = zeros(d{1}(2),d{1}(1));
        imshow(I,'border', 'tight');
        pb = pb/max(max(pb));
        draw_bry_contours_prune_pb(d{2}, pb);
        h = figure(1);
    %     print_pdf([src(1:end-3), 'pdf'], h);
        set(h, 'PaperPositionMode','auto')
        print(h,'-dpng','-r0',[src(1:end-4),'.png']);
%         movefile([src(1:end-4),'.png'], [src(1:end-6),'.png'])
    end

% end