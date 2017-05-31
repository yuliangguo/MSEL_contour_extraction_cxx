clear all; close all;

src_path = '/home/guoy/Desktop/visualization/amsterdam_house_contours_orig_size/PB_Kokkinos/';

cem_files = dir([src_path '*.cem']);

for c = 1:length(cem_files)
    c
    filename = [src_path cem_files(c).name];
    [cem, edg] = load_contours(filename);
    edgename = [filename(1:end-4) '.edg'];
    save_edg(edgename, edg, cem{1});
end