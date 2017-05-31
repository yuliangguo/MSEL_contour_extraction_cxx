src_path = '/home/guoy/Desktop/BSR_yuliang/CFGD/coarse/s3/';
src_files = dir([src_path, '*.cem']);

for i = 1:length(src_files)

file_name = [src_path, src_files(i).name];   
d2=load_contours(file_name);
I = zeros(d2{1}(2),d2{1}(1));

close all;
h = figure(1);
set(h, 'PaperPositionMode','auto')
imshow(I, 'border', 'tight');
hold on;

draw_contours(d2{2},0,0,[0.99,0.99,0.99]) % draw binary map

outputname = [file_name(1:end-4), '.png'];
print(h,'-dpng', '-r0', outputname); % never use .jpg, -r0 to keep the resolution

end