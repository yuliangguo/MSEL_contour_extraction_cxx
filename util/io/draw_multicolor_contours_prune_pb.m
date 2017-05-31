function draw_multicolor_contours_prune_pb(contours, pb_map, pb_th)

 if (nargin<3), pb_th = 0; end
% if (nargin<3), rand=1; end
% if (nargin<4), col = [0 0 0]; end

con_cnt = length(contours);

% this is for multicolor
colourmp = hsv(con_cnt);    % HSV colour map with con_cnt entries
colourmp = colourmp(randperm(con_cnt),:);  % Random permutation
% 
% % this is for white
% colourmp(:,1) = 0.99;
% colourmp(:,2) = 0.99;
% colourmp(:,3) = 0.99;

% % this is for red
% colourmp(:,1) = 255/255;
% colourmp(:,2) = 0/255;
% colourmp(:,3) = 0/255;

for i = 1:con_cnt
    c_pb = contour_pb(contours{1,i}, pb_map); % pb of contour is the avg value of pb along its edges
    
    if(c_pb>pb_th)
%         col = ones(1,3)*c_pb;
        line(contours{i}(:,1)+1, contours{i}(:,2)+1,'color', colourmp(i,:), 'LineWidth', 2);
    end

end