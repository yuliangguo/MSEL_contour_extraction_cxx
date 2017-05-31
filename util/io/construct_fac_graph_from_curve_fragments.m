function G = construct_fac_graph_from_curve_fragments (edges, cf_idx, contours)
% construct a graph G: graph node represents connecting curve fragments end
% points, graph factor represents each curve fragments

G = init_graph();

used_edge_ids = [];
num_of_short_cfs = 0;
for c = 1:length(cf_idx)
    
%     % ignore the extremely short curve fragments
%     if(length(cf_idx{c})<=2)
%         num_of_short_cfs = num_of_short_cfs+1;
% %         continue;        
%     end
    
    
    search_1 = find(used_edge_ids==cf_idx{c}(1));
    if(isempty(search_1))
        [G, id_1] = add_varnode(G, cf_idx{c}(1));
        used_edge_ids = [used_edge_ids cf_idx{c}(1)];
    else
        id_1 = search_1(1);
    end
    
    search_2 = find(used_edge_ids==cf_idx{c}(end));
    if(isempty(search_2))
        [G, id_2]= add_varnode(G, cf_idx{c}(end));
        used_edge_ids = [used_edge_ids cf_idx{c}(end)];
    else
        id_2 = search_2(1);
    end
    
    G = add_facnode(G, contours{c}, [id_1 id_2]);
end

% num_of_short_cfs

end