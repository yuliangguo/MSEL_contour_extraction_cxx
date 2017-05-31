function draw_contours_sorted_col(contours, col)

con_cnt = length(contours);

for i = 1:con_cnt
    line(contours{i}(:,1)+1, contours{i}(:,2)+1,'color', col(i,:), 'LineWidth', 1);
    
    plot(contours{i}(1,1)+1, contours{i}(1,2)+1, 'y.', 'MarkerSize', 3);
    plot(contours{i}(end,1)+1, contours{i}(end,2)+1, 'y.', 'MarkerSize', 3);
end
