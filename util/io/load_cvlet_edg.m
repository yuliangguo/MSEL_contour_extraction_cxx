function edg = load_cvlet_edg(filename)

fid = fopen(filename);

edge_cnt = 1;
width = 0;
height=0;

while 1
    lineBuffer = fgetl(fid);
    if ~ischar(lineBuffer), break, end

    %ignore other comment lines and empty lines
    if (length(lineBuffer)<2 | lineBuffer(1)=='#')
        continue;
    end

    if (strncmp(lineBuffer, 'WIDTH=', length('WIDTH=')))
        width = strread(lineBuffer,'WIDTH=%d');
        continue;
    elseif (strncmp(lineBuffer, ' WIDTH=', length(' WIDTH=')))
        width = strread(lineBuffer,' WIDTH=%d');
        continue;
    end
    
    if (strncmp(lineBuffer, 'HEIGHT=', length('HEIGHT=')))
        height = strread(lineBuffer,'HEIGHT=%d');
        continue;
    elseif (strncmp(lineBuffer, ' HEIGHT=', length(' HEIGHT=')))
        height = strread(lineBuffer,' HEIGHT=%d');
        continue;
    end
    
    % read the line with the edge count info
    if (strncmp(lineBuffer, 'EDGE_COUNT=', length('EDGE_COUNT=')))
        num_edges = strread(lineBuffer,'EDGE_COUNT=%d');
            num_edges = strread(lineBuffer,'EDGE_COUNT=%d');
            edg = zeros(num_edges, 5);
        continue;
    end

    if(strncmp(lineBuffer, '[BEGIN EDGEMAP]', length('[BEGIN EDGEMAP]')))
        lineBuffer = fgetl(fid);
        while 1
            lineBuffer = fgetl(fid);
            if(strncmp(lineBuffer, '[END EDGEMAP]', length('[END EDGEMAP]')))
               break; 
            end
            [unused, x, y, dir, conf] = strread(lineBuffer,' [%d] [%f, %f]   %f %f');
            edg(edge_cnt,:) = [x y dir conf 0];
            edge_cnt = edge_cnt + 1; 
        end
    else
        continue;
    end
    
    if(strncmp(lineBuffer, '[BEGIN CVLETMAP]', length('[BEGIN CVLETMAP]')))
        break;
    end
end
    %close the file
    fclose(fid);

    % update edge coordinates to matlab coordinates
    edg(:,1:2) = edg(:,1:2)+1;
end
