function write_cem(filename, contours, h, w)


    edgelist = contours;
    edge_id_list = cell(size(edgelist));

    edge_map = edgelist{1,1};
    count = size(edgelist{1,1}, 1);
    edge_id_list{1,1} = 0:(count -1);
    for i = 2:size(edgelist,2)
        length = size(edgelist{1,i}, 1);
        edge_map = [edge_map; edgelist{1,i}];
        edge_id_list{1,i} = count:(count + length -1);
        count = count + length;
    end

% write new contour map files
    fprintf('saving into .cem file\n');
    fout = fopen(filename, 'w');
    fprintf(fout, '.CEM v2.0\n');
    fprintf(fout, ['size=[',num2str(w),' ', num2str(h), ']\n']);
    fprintf(fout, '[Edgemap]\n');
    fprintf(fout, 'count=%d\n', size(edge_map,1));
    for j=1:size(edge_map,1)
      % when saving edges, we make posions -1, because in cxx, index start
      % from 0, but in matlab index starts from 1
      s = ['(',num2str(edge_map(j,1)),', ', num2str(edge_map(j,2)),')  ',num2str(edge_map(j,3)),'   ',num2str(edge_map(j,4)),'   0'];
      fprintf(fout, '%s\n', s);
    end
    fprintf(fout, '[Contours]\n');
    fprintf(fout, 'count=%d\n', size(edgelist,2));
    for j = 1:size(edge_id_list,2)
        contour = edge_id_list{1,j};
        contour_str = num2str(contour(1,1));
        for k = 2:size(contour, 2)
            contour_str = [contour_str,' ', num2str(contour(1,k))];
        end
        fprintf(fout, '[%s ]\n', contour_str);
    end
    fprintf(fout, '[Contour Properties]\n');
    fprintf(fout, '# <len> <avg. str> <mean con> <Lstd> <Rstd> <avg. d2f> <avg. k> <max k>\n');
    for j = 1:size(edge_id_list,2)
      s = [num2str(size(edge_id_list{1,j},2)),' 2 0 0 0 0 0.5 0.5'];
      fprintf(fout, '%s\n', s);
    end
    
    fclose(fout);