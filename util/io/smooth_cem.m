function smooth_cem(src_file, dst_file)

    psi = 1; %this is the step length, default as 1;
    ns = 1; %iterations of smmothing, default as 1;

    d2=load_contours_gt_smooth(src_file);

    % smooth each contour and replace with new positions
    % compute new orientation of each edge, based on new positions
    fprintf('smoothing contours\n');
    for i = 1:size(d2{2,1},2)
        edges = d2{2,1}{1,i};
%         if(size(edges,1)<3)
%             continue;
%         end
        new_edges = edges;
        n_edgs = size(edges,1);
        for k = 1:ns
            for j = 2:n_edgs-1
                middle_p = edges(j,1:2);
                left_p = edges(j-1,1:2);
                right_p = edges(j+1,1:2);

                v1 = left_p - middle_p;
                v2 = right_p - middle_p;

                nv1 = norm(v1);
                nv2 = norm(v2);

                v1n = v1/nv1;
                v2n = v2/nv2;

                if(nv1 < nv2)
                    nrm = (v1 + v2n*nv1)/4.0;
                else
                    nrm = (v2 + v1n*nv2)/4.0;
                end

                cs_i = middle_p + psi * nrm;
                new_edges(j,1:2) = cs_i;
            end
            edges = new_edges;
        end

        % recompute the orientations for each edge within this contours
        for j = 2:n_edgs-1
            left_p = edges(j-1,1:2);
            right_p = edges(j+1,1:2);

            delta = right_p - left_p;
            edges(j,3) = atan2(delta(1,2), delta(1,1));
            if(edges(j,3)<0)
                %change the interval from [-pi, pi) to [0,2*pi)
                edges(j,3) = edges(j,3) + 2*pi;
            end
        end

        % deal with the orientations on both ends
        edges(1,3) = edges(2,3);
        edges(n_edgs,3) = edges(n_edgs - 1,3);

        % assign new contour back into contour map
        d2{2,1}{1,i} = edges;

    end

    %put all edges into one edgemap, and build up a contour list with edge
    %index
    edgelist = d2{2,1};
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
%     dst = [dst_path, src_files(f).name(1:end-4) '.cem'];
    fout = fopen(dst_file, 'w');
    fprintf(fout, '.CEM v2.0\n');
    fprintf(fout, ['size=[',num2str(d2{1,1}(1,1)),' ', num2str(d2{1,1}(1,2)), ']\n']);
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