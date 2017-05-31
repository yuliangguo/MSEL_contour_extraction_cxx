clear all; close all;
src_path = '/home/guoy/Desktop/visualization/amsterdam_house_contours_orig_size/PB_Kokkinos/';
dst_path = '/home/guoy/Desktop/visualization/amsterdam_house_contours_orig_size/PB_Kokkinos/';

cem_files = dir([src_path, '*.cem']);

for i = 1:length(cem_files)
    i
    input_file = [src_path cem_files(i).name];
    CEM = load_contours(input_file);
    det_save_cemv([input_file(1:end-4), '.cemv'], CEM{2});
%     keyboard;
end