clear all; close all;
img_path = '/home/guoy/Desktop/CFGD/GT_img/';

images = dir([img_path '*jpg']);
for i = 1:length(images)
    img_name = images(i).name;
    src = [img_path, img_name];
    img = imread(src);        
    resize_img = imresize(img, 0.5);
    
%     dst = [img_path, 'half_scale/', img_name];
%     imwrite(resize_img, dst, 'JPEG');
    
    dst = [img_path, 'half_scale/', img_name(1:end-3), 'png'];
    imwrite(resize_img, dst, 'PNG');
end
