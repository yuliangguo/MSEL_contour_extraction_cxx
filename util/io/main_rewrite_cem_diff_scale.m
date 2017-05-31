clear all; close all;
src_path = '/home/guoy/Desktop/visualization/amsterdam_house_contours/Kovesi/';
dst_path = '/home/guoy/Desktop/visualization/amsterdam_house_contours_orig_size/PB_Kovesi/';
mkdir(dst_path)
files = dir([src_path '*cem']);
scale = 4;

for f = 1:length(files)
    f
    filename = [src_path files(f).name];
    dst = [dst_path files(f).name];
    tic;
    fid = fopen(filename);
    fout = fopen(dst, 'w');

    if (fid<0)
        disp('File not found.');
        return;
    end

    lineBuffer = fgetl(fid); %read in the first line

    if (~strncmp(lineBuffer, '.CEM v2.0', length('.CEM v2.0')))
        return; %wrong version
    end


    %scan the file to read the contour information
    while 1,
        lineBuffer = fgetl(fid);
        if ~ischar(lineBuffer), break, end

        %ignore comment lines and empty lines
        if (length(lineBuffer)<2)
            continue;
        end

        % read the line with the edgemap size info
        if (strncmp(lineBuffer, 'size=', length('size=')))
            [w, h] = strread(lineBuffer,'size=[%d %d]');
            fprintf(fout, '.CEM v2.0\n');
            fprintf(fout, ['size=[',num2str(w),' ',num2str(h), ']\n']);
            continue;
        end

        %read the edgemap block
        if (strncmp(lineBuffer, '[Edgemap]', length('[Edgemap]')))
            %read the next line with the edge count
            lineBuffer = fgetl(fid);
            num_edges = strread(lineBuffer,'count=%d');
            fprintf(fout, '[Edgemap]\n');
            fprintf(fout, 'count=%d\n', num_edges);
            %read in all the edges
            for i=1:num_edges,
                lineBuffer = fgetl(fid);
                [x, y, dir, conf, d2f] = strread(lineBuffer,'(%f, %f)\t%f\t%f\t%f');
    %             edg(i,:) = [x, y, dir, conf, d2f];
                s = ['(',num2str(x*scale),', ', num2str(y*scale),')  ',num2str(dir),' ', num2str(conf), ' ', num2str(d2f)];
                fprintf(fout, '%s\n', s);           
            end
            continue;
        end

        % read the contours block
        if (strncmp(lineBuffer, '[Contours]', length('[Contours]')))
            %read the next line with the contour count
            lineBuffer = fgetl(fid);
            num_contours = strread(lineBuffer,'count=%d');
            fprintf(fout, '[Contours]\n');
            fprintf(fout, 'count=%d\n', num_contours);

            %read in all the contours
            for i=1:num_contours,
                %read in the list of edge ids
                lineBuffer = fgetl(fid);
                fprintf(fout, '%s\n', lineBuffer);

    %             e_ids = strread(lineBuffer(2:length(lineBuffer)-1),'%d','delimiter',' ')';
    %   
    %             chain = [];
    %             for j=1:length(e_ids)
    %                 chain(j,:) = edg(e_ids(j)+1,:);
    %             end
    % 
    %             %store chain
    %             contours{i} = chain;
            end

    %         %save the contours
    %         cem{2} = contours;
            continue;
        end

        % read the contour properties block
    %     if (strncmp(lineBuffer, '[Contour Properties]', length('[Contour Properties]')))
    %         lineBuffer = fgetl(fid);
            fprintf(fout, '%s\n', lineBuffer);

    %         %read in all the contour properties
    %         con_props = [];
    %         for i=1:num_contours,
    %             %read in the list of contour properties
    %             lineBuffer = fgetl(fid);
    %             
    %             %<len> <avg. str> <mean con> <Lstd> <Rstd> <avg. d2f> <avg. k> <max k>
    %             con_props(i,:) = strread(lineBuffer,'%f','delimiter',' ');
    %         end
    %         
    %         %save the properties
    %         cem{3} = con_props;
    %     end
    end

    %close the file
    fclose(fid);
    fclose(fout);
    toc;
end