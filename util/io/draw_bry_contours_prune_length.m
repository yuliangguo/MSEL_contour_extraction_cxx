function draw_bry_contours_prune_length(contours, length_thresh, rand, col)

if (nargin<2), length_thresh = 0; end
if (nargin<3), rand=1; end
if (nargin<4), col = [0 0 0]; end

con_cnt = length(contours);

% this is for multicolor
colourmp = hsv(con_cnt);    % HSV colour map with con_cnt entries
colourmp = colourmp(randperm(con_cnt),:);  % Random permutation

% this is for white
colourmp(:,1) = 0.99;
colourmp(:,2) = 0.99;
colourmp(:,3) = 0.99;

% % this is for red
% colourmp(:,1) = 255/255;
% colourmp(:,2) = 0/255;
% colourmp(:,3) = 0/255;

for i = 1:con_cnt
    len = contour_length_mex(contours{1,i}');
    if (len<length_thresh)
        continue;
    end

    if (rand==1)
        line(contours{i}(:,1)+1, contours{i}(:,2)+1,'color',colourmp(i,:), 'LineWidth', 2);
    else

        line(contours{i}(:,1)+1, contours{i}(:,2)+1,'color', col);

    end
    clear len;
end