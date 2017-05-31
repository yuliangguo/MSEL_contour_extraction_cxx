function draw_contours(contours, thresh, rand, col)

if (nargin<2), thresh = 0; end
if (nargin<3), rand=1; end
if (nargin<4), col = [0 0 0]; end

con_cnt = length(contours);

% this is for multicolor
colourmp = hsv(con_cnt);    % HSV colour map with con_cnt entries
colourmp = colourmp(randperm(con_cnt),:);  % Random permutation

% this is for green
% colourmp(:,1) = 0/255;
% colourmp(:,2) = 255/255;
% colourmp(:,3) = 0/255;

% % this is for red
% colourmp(:,1) = 255/255;
% colourmp(:,2) = 0/255;
% colourmp(:,3) = 0/255;

for i = 1:con_cnt
    if (size(contours{i},1)<thresh)
        continue;
    end

    if (rand==1)
        line(contours{i}(:,1)+1, contours{i}(:,2)+1,'color',colourmp(i,:), 'LineWidth', 2);
    else

        line(contours{i}(:,1)+1, contours{i}(:,2)+1,'color', col, 'LineWidth', 2);

    end
    
    plot(contours{i}(1,1)+1, contours{i}(1,2)+1, 'y.', 'MarkerSize', 4);
    plot(contours{i}(end,1)+1, contours{i}(end,2)+1, 'y.', 'MarkerSize', 4);
end
