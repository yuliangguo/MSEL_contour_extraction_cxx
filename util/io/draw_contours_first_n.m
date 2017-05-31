function draw_contours_first_n(contours, n, thresh, rand, col)

if (nargin<3), thresh = 0; end
if (nargin<4), rand=1; end
if (nargin<5), col = [0 0 0]; end

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

for i = 1:min(n,con_cnt)
    if (size(contours{i},1)<thresh)
        continue;
    end

    if (rand==1)
        line(contours{i}(:,1)+1, contours{i}(:,2)+1,'color',colourmp(i,:), 'LineWidth', 2);
    else

        line(contours{i}(:,1)+1, contours{i}(:,2)+1,'color', col, 'LineWidth', 1);

    end
end
