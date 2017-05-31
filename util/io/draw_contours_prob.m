function draw_contours_prob(contours, prob_vec, thresh)

if (nargin<3), thresh = 0; end


con_cnt = length(contours);

for i = 1:con_cnt
    if (size(contours{i},1)<thresh)
        continue;
    end

    col = [1,1,1]*prob_vec(i);

    line(contours{i}(:,1)+1, contours{i}(:,2)+1,'color', col, 'LineWidth', 1);
    
%     plot(contours{i}(1,1)+1, contours{i}(1,2)+1, 'y.', 'MarkerSize', 3);
%     plot(contours{i}(end,1)+1, contours{i}(end,2)+1, 'y.', 'MarkerSize', 3);
end
