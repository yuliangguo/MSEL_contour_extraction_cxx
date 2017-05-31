function c_pb = contour_pb(contour, pb_map)

pb_list = [];

for e = 1: size(contour,1)
    pb = pb_map(round(contour(e,2)+1), round(contour(e,1)+1));
    pb_list = [pb_list pb];
end

c_pb = mean(pb_list);

end