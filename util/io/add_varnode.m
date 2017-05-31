function [G, id] = add_varnode( G, actual_edge_id )
% ADD_VARNODE - Add variable node to factor graph 'G'.
%
  id = numel(G.var) + 1;
  
%   v.name = name;
%   v.dim = numel(vals);
%   v.vals = vals;
  v.dim = [];
  v.id = id;
  v.nbrs_fac = [];
  v.p = [];
  v.merged = 0;
%   v.observed = 0;
  v.gt_label = -1; %% -1 means undetermined, 0: break, 1:merge
  v.actual_edge_id = actual_edge_id;
  G.var = cat( 1, G.var, v );

end