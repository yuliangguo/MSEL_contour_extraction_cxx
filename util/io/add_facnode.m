function G = add_facnode(G, contour, nbrs_var)
% ADD_FACNODE - Add factor node to factor graph 'G'.  
%
% INPUTS:
%   G - Factor graph
%
%   p - Potential matrix with p(i,j,...) the potential for
%       x_a=i, x_b=j, ...
%
%   varargin - Variable nodes involved in factor.  Order matches dimensions
%              potential, e.g. varargin = { a, b, ... } => p(x_a,x_b,...).
%
%   f.p = p;

  f.nbrs_var = nbrs_var;
  f.id = numel(G.fac) + 1;
%   f.shock_edge_id = shock_edge_id;
  f.contour = contour; % save the curve fragment
  f.removed = 0; % =1 when some of the curve fragment is merged into another or pruned 
%   % check dimensions
%   if ( ( numel(f.nbrs_var) > 1 ) && ( numel(f.nbrs_var) ~= ndims(p)) ) || ...
%     ( ( numel(f.nbrs_var) == 1 ) && ( ndims(p) ~= 2 ) )
%     error('add_facnode: Factor dimensions does not match size of domain.');
%   end
  
  G.fac = cat( 1, G.fac, f );
  for I=f.nbrs_var
    G.var(I).nbrs_fac = [ G.var(I).nbrs_fac; f.id ];
    G.var(I).dim = length(G.var(I).nbrs_fac);
  end
end