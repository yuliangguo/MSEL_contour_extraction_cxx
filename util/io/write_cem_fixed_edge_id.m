function write_cem_fixed_edge_id(filename, contours, edges, h, w, cfrags_idx)


    edgelist = contours;
    edge_id_list = cell(size(edgelist));

    edge_map = edges;
    
    if(nargin>=6)  % if the curve frags with idx are given, use it directly
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
        fprintf(fout, 'count=%d\n', length(cfrags_idx));
        for j = 1:length(cfrags_idx)
            contour = cfrags_idx{j}-1; %%%%%%%%%%%%%%%%%%%  change the id to cxx coord
            contour_str = mat2str(contour);
            fprintf(fout, '%s\n', contour_str);
        end
        fprintf(fout, '[Contour Properties]\n');
        fprintf(fout, '# <len> <avg. str> <mean con> <Lstd> <Rstd> <avg. d2f> <avg. k> <max k>\n');
        for j = 1:size(edge_id_list,2)
          s = [num2str(size(edge_id_list{1,j},2)),' 2 0 0 0 0 0.5 0.5'];
          fprintf(fout, '%s\n', s);
        end

        fclose(fout);
        
    else  %  construct hashtable to save the cfrags with edges
    
    
        values = cell(1,size(edges,1));
        Keys = cell(1,size(edges,1));

        for i = 1:size(edges,1)
            Key = mat2str(edges(i,1:3));
            values{i} = i-1; %%%%%%%%%%%%%%%%%%%  change the id to cxx coord
            Keys{i} = Key;
        end

        Hashtable = containers.Map(Keys, values); 

    %     count = size(edges, 1);
        for i = 1:size(edgelist,2)
            c_edges = edgelist{1,i};
            length_c = size(c_edges, 1);
            edge_id_list{1,i} = zeros(1,length_c);
            for j = 1:length_c
                edge_id_list{1,i}(j) = Hashtable(mat2str(c_edges(j,1:3)));
            end
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
    
    end