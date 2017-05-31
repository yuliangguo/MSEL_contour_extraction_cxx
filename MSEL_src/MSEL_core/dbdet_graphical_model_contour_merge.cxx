#include <vcl_map.h>
#include "dbgl_diffgeom.h"
#include "dbdet_graphical_model_contour_merge.h"
#include <vcl_algorithm.h>

struct var_node {
  /*enum label { UNDETERMINED, BREAK, MERGE };
  unsigned id;*/
  unsigned dim;
  unsigned edgel_id;
  vcl_vector<unsigned> n_facs;
  /*bool merged;
  double p;
  label gt_label;*/

  var_node(/*unsigned uid,*/ unsigned e_id) : /*id(uid),*/ dim(0), edgel_id(e_id)/*, merged(false), p(0.0), gt_label(UNDETERMINED)*/ {} 

  void push_fac(unsigned fac_id)
  {
    n_facs.push_back(fac_id);
    dim = n_facs.size();
  }
};

struct fac_node {
  dbdet_edgel_chain * chain;
  //unsigned id;
  vcl_vector<unsigned> n_vars;
  //bool removed;

  fac_node(dbdet_edgel_chain * c, /*unsigned uid,*/ vcl_vector<unsigned> vars) : chain(c), n_vars(vars)/*, id(uid), removed(false) */{}
};

class dbdet_graphical_model_contour_merge::dbdet_factor_graph {
  
public:
  vcl_vector<var_node> var;
  vcl_vector<fac_node> fac;
 
  dbdet_factor_graph(dbdet_curve_fragment_graph & CFG)
  {
    var.reserve(2 * CFG.frags.size());
    fac.reserve(CFG.frags.size());
    typedef vcl_pair<vcl_map<unsigned, unsigned>::iterator, bool> insert_ret;
    vcl_map<unsigned, unsigned> used_edgels;// edgel_id, var_node_id

    for (dbdet_edgel_chain_list_const_iter it=CFG.frags.begin(); it != CFG.frags.end(); it++)
    {
      dbdet_edgel * e1 = (*it)->edgels.front();
      dbdet_edgel * e2 = (*it)->edgels.back();

      vcl_vector<unsigned> var_ids;
      unsigned id = var.size();
      insert_ret ret = used_edgels.insert(vcl_pair<unsigned, unsigned>(e1->id, id));
      if(ret.second)
      {
        add_var_node(e1->id);
        var_ids.push_back(id);
      }
      else
      {
        var_ids.push_back((ret.first)->second);
      }

      id = var.size();
      ret = used_edgels.insert(vcl_pair<unsigned, unsigned>(e2->id, id));
      if(ret.second)
      {
        add_var_node(e2->id);
        var_ids.push_back(id);
      }
      else
      {
        var_ids.push_back((ret.first)->second);
      }
      add_fac_node(*it, var_ids);
    }
  }

private:
  void add_var_node(unsigned edgel_id)
  {
    var.push_back(var_node(/*var.size() + 1, */edgel_id));
  }

  void add_fac_node(dbdet_edgel_chain * c, vcl_vector<unsigned> vars)
  {
    unsigned id = fac.size();
    fac.push_back(fac_node(c,/* id,*/ vars));
    for (unsigned i = 0; i < vars.size(); ++i)
    {
       var[vars[i]].push_fac(id);
    }
  }
};

double compute_merge_prob_sem(
      y_feature_vector cues,
      const y_params_0_vector & beta0,
      const y_params_0_vector & fmean0
      )
{
  double sum = 0.0; 
  for (unsigned i = 0; i < cues.size(); ++i)
    sum += (cues[i] - fmean0[i]) * beta0[i];
  return 1.0 / (1.0 + exp(-sum));
}

void dbdet_graphical_model_contour_merge::
dbdet_merge_contour(
      dbdet_curve_fragment_graph & CFG,
      y_params_1_vector & beta1,
      y_params_1_vector & fmean1,
      y_params_0_vector & beta0,
      y_params_0_vector & fmean0,
      dbdet_curve_fragment_graph & newCFG
      )
{
  deep_copy_cfg(CFG, newCFG);
  dbdet_factor_graph g(newCFG);

  for (unsigned k = 0; k < g.var.size(); ++k)
  {
    var_node & cur_node = g.var[k];

    if (cur_node.dim == 2)
    {
      unsigned c1_id = cur_node.n_facs.front();
      unsigned c2_id = cur_node.n_facs.back();

      if (c1_id == c2_id)
        continue;

      dbdet_edgel_chain * c1 = g.fac[c1_id].chain;
      dbdet_edgel_chain * c2 = g.fac[c2_id].chain;

      dbdet_edgel_chain c1_cut, c2_cut;
      
      //Original direction c1 -> node -> c2
      //Our directions c1 -> node <- c2

      unsigned cut_size = vcl_min(nbr_num_edges, static_cast<unsigned>(c1->edgels.size()));
      c1_cut.edgels.resize(cut_size);
      if (c1->edgels.front()->id == cur_node.edgel_id)
      {
        for (unsigned i = 0; i < cut_size; ++i)
          c1_cut.edgels[i] = c1->edgels[cut_size - i - 1];
      }
      else /*if(c1->edgels.back()->id == cur_node.edgel_id)*/
      {
        unsigned start_i = c1->edgels.size() - cut_size;
        for (unsigned i = 0; i < cut_size; ++i)
          c1_cut.edgels[i] = c1->edgels[start_i + i];
      }
      cut_size = vcl_min(nbr_num_edges, static_cast<unsigned>(c2->edgels.size()));
      c2_cut.edgels.resize(cut_size);

      if (c2->edgels.front()->id == cur_node.edgel_id)
      {
        for (unsigned i = 0; i < cut_size; ++i)
          c2_cut.edgels[i] = c2->edgels[cut_size - i - 1];
      }
      else /*if(c2->edgels.back()->id == cur_node.edgel_id)*/
      {
        unsigned start_i = c2->edgels.size() - cut_size;
        for (unsigned i = 0; i < cut_size; ++i)
          c2_cut.edgels[i] = c2->edgels[start_i + i];
      }
      y_feature_vector c1_features, c2_features;
      cues_computer.compute_all_cues(c1_cut, &c1_features);
      cues_computer.compute_all_cues(c2_cut, &c2_features);

      double c1_len = c1_features[y_features::Y_LEN];
      double c2_len = c2_features[y_features::Y_LEN];

      y_params_0_vector cues;
      cues[y_features::Y_ONE] = 1.0;

      //Change sign because of inverted (converging) direction
      c2_features[y_features::Y_BG_GRAD] *= -1.;
      c2_features[y_features::Y_SAT_GRAD] *= -1.;
      c2_features[y_features::Y_HUE_GRAD] *= -1.;

      for (unsigned i = 1; i < cues.size(); ++i)
      {
        cues[i] = vcl_abs(c1_features[i] - c2_features[i]);
      }

      double geom_diff, texture_diff;

      dbdet_degree_2_node_cues(c1_cut, c2_cut, geom_diff, texture_diff);
      cues[y_params_0::Y_GEOM] = geom_diff;
      cues[y_params_0::Y_TEXTURE] = texture_diff;

      bool merge = false;

      //for very short curve, just decide merging based on geometry;
      if (c1_len < dbdet_yuliang_const::nbr_len_th || c2_len < dbdet_yuliang_const::nbr_len_th)
      {
        double p = 1.0 / (1.0 + exp(-((1.0 - fmean1[0]) * beta1[0] + (cues[y_params_0::Y_GEOM] - fmean1[1]) * beta1[1])));
        merge = p > dbdet_yuliang_const::merge_th_geom ? true : false;
      }
      else
      {
        double p = compute_merge_prob_sem(cues, beta0, fmean0);
        merge = p > dbdet_yuliang_const::merge_th_sem ? true : false;
      }

      if (merge)
      {
        dbdet_merge_at_degree_2_node(g, c1_id, c2_id, k, cur_node.edgel_id);
      }
    }
  }

  //dim 3 node should be dealt with after all dim 2 node are solved 
  for (unsigned k = 0; k < g.var.size(); ++k)
  {
    var_node & cur_node = g.var[k];
    if (cur_node.dim == 3)
    {
      unsigned c1_id = cur_node.n_facs[0];
      unsigned c2_id = cur_node.n_facs[1];
      unsigned c3_id = cur_node.n_facs[2];
      if (c1_id == c2_id || c1_id == c3_id || c2_id == c3_id)
        continue;

      dbdet_edgel_chain * c1 = g.fac[c1_id].chain;
      dbdet_edgel_chain * c2 = g.fac[c2_id].chain;
      dbdet_edgel_chain * c3 = g.fac[c3_id].chain;

      dbdet_edgel_chain c1_cut, c2_cut, c3_cut;
      // always make directions c1 -> node <- c2
      //                               ^
      //                               |
      //                               c3

      unsigned cut_size = vcl_min(nbr_num_edges, static_cast<unsigned>(c1->edgels.size()));
      c1_cut.edgels.resize(cut_size);
      if (c1->edgels.front()->id == cur_node.edgel_id)
      {
        for (unsigned i = 0; i < cut_size; ++i)
          c1_cut.edgels[i] = c1->edgels[cut_size - i - 1];
      }
      else /*if(c1->edgels.back()->id == cur_node.edgel_id)*/
      {
        unsigned start_i = c1->edgels.size() - cut_size;
        for (unsigned i = 0; i < cut_size; ++i)
          c1_cut.edgels[i] = c1->edgels[start_i + i];
      }

      cut_size = vcl_min(nbr_num_edges, static_cast<unsigned>(c2->edgels.size()));
      c2_cut.edgels.resize(cut_size);
      if (c2->edgels.front()->id == cur_node.edgel_id)
      {
        for (unsigned i = 0; i < cut_size; ++i)
          c2_cut.edgels[i] = c2->edgels[cut_size - i - 1];
      }
      else /*if(c2->edgels.back()->id == cur_node.edgel_id)*/
      {
        unsigned start_i = c2->edgels.size() - cut_size;
        for (unsigned i = 0; i < cut_size; ++i)
          c2_cut.edgels[i] = c2->edgels[start_i + i];
      }

      cut_size = vcl_min(nbr_num_edges, static_cast<unsigned>(c3->edgels.size()));
      c3_cut.edgels.resize(cut_size);
      if (c3->edgels.front()->id == cur_node.edgel_id)
      {
        for (unsigned i = 0; i < cut_size; ++i)
          c3_cut.edgels[i] = c3->edgels[cut_size - i - 1];
      }
      else /*if(c3->edgels.back()->id == cur_node.edgel_id)*/
      {
        unsigned start_i = c3->edgels.size() - cut_size;
        for (unsigned i = 0; i < cut_size; ++i)
          c3_cut.edgels[i] = c3->edgels[start_i + i];
      }

      y_feature_vector c1_features, c2_features, c3_features;
      cues_computer.compute_all_cues(c1_cut, &c1_features);
      cues_computer.compute_all_cues(c2_cut, &c2_features);
      cues_computer.compute_all_cues(c3_cut, &c3_features);

      y_params_0_vector cues_12, cues_13, cues_23;
      cues_12[y_params_0::Y_ONE] = cues_13[y_params_0::Y_ONE] = cues_23[y_params_0::Y_ONE] = 1.0;

      for (unsigned i = y_params_0::Y_ABS_K; i < cues_12.size(); ++i)
      {
        cues_12[i] = vcl_abs(c1_features[i] - c2_features[i]);
        cues_13[i] = vcl_abs(c1_features[i] - c3_features[i]);
        cues_23[i] = vcl_abs(c2_features[i] - c3_features[i]);
      }

      //Sums because of inverted (converging) direction
      for (unsigned i = y_params_0::Y_BG_GRAD; i <= y_params_0::Y_HUE_GRAD; ++i)
      {
        cues_12[i] = vcl_abs(c1_features[i] + c2_features[i]);
        cues_13[i] = vcl_abs(c1_features[i] + c3_features[i]);
        cues_23[i] = vcl_abs(c2_features[i] + c3_features[i]);
      }


      dbdet_degree_2_node_cues(c1_cut, c2_cut, cues_12[y_params_0::Y_GEOM], cues_12[y_params_0::Y_TEXTURE]);
      dbdet_degree_2_node_cues(c1_cut, c3_cut, cues_13[y_params_0::Y_GEOM], cues_13[y_params_0::Y_TEXTURE]);
      dbdet_degree_2_node_cues(c2_cut, c3_cut, cues_23[y_params_0::Y_GEOM], cues_23[y_params_0::Y_TEXTURE]);

      double p_12 = compute_merge_prob_sem(cues_12, beta0, fmean0);
      double p_13 = compute_merge_prob_sem(cues_13, beta0, fmean0);
      double p_23 = compute_merge_prob_sem(cues_23, beta0, fmean0);

      double max = p_12 > p_13 ? (p_12 > p_23 ? p_12 : p_23) : (p_13 > p_23 ? p_13 : p_23);

      if(max > dbdet_yuliang_const::merge_th_sem)
      {
        if (max == p_12)
        {
          dbdet_merge_at_degree_2_node(g, c1_id, c2_id, k, cur_node.edgel_id);
        }
        else if (max == p_13)
        {
          dbdet_merge_at_degree_2_node(g, c1_id, c3_id, k, cur_node.edgel_id);
        }
        else
        {
          dbdet_merge_at_degree_2_node(g, c2_id, c3_id, k, cur_node.edgel_id);
        }
      }
    }
  }

  for (dbdet_edgel_chain_list_iter it=newCFG.frags.begin(); it != newCFG.frags.end();)
  {
    if((*it)->edgels.size() == 0)
      newCFG.frags.erase(it++);
    else
      it++;
  }
}

void dbdet_graphical_model_contour_merge::
dbdet_degree_2_node_cues(
      dbdet_edgel_chain & c1,
      dbdet_edgel_chain & c2,
      double & geom_diff,
      double & tex_diff,
      bool invert
      )
{
  unsigned nbr_range_th = vcl_min(c1.edgels.size(), c2.edgels.size());
  vgl_vector_2d<double> a_ori = c1.edgels.back()->pt - c1.edgels[c1.edgels.size() - nbr_range_th]->pt;
  vgl_vector_2d<double> b_ori = c2.edgels.back()->pt - c2.edgels[c2.edgels.size() - nbr_range_th]->pt;
  
  geom_diff = (b_ori.x() * a_ori.x() + b_ori.y() * a_ori.y()) / 
        (vcl_sqrt(b_ori.x() * b_ori.x() + b_ori.y() * b_ori.y()) * vcl_sqrt(a_ori.x() * a_ori.x() + a_ori.y() * a_ori.y()));

  if(invert)
    geom_diff *= -1.0;

  y_hist_vector left_1, left_2, right_1, right_2;
  dbgl_compute_normals(c1, &temp_n);
  compute_texture_hist(c1, temp_n, left_1, right_1);

  dbgl_compute_normals(c2, &temp_n);

  if (invert)
    compute_texture_hist(c2, temp_n, right_2, left_2);
  else
    compute_texture_hist(c2, temp_n, left_2, right_2);

  tex_diff = 0.0;

  //chisq dist
  for (unsigned k = 0; k < y_hist_size; ++k)
  {
    double v1 = left_1[k];
    double v2 = left_2[k];
    double den = ((v1 + v2) == 0.0) ? 1e-16 : v1 + v2;
    tex_diff += ((v1 - v2) * (v1 - v2)) / den;

    v1 = right_1[k];
    v2 = right_2[k];
    den = ((v1 + v2) == 0.0) ? 1e-16 : v1 + v2;
    tex_diff += ((v1 - v2) * (v1 - v2)) / den;
  }
  tex_diff /= 2.0;
}

void dbdet_graphical_model_contour_merge::
dbdet_merge_at_degree_2_node(
      dbdet_factor_graph & G,
      unsigned c1_id,
      unsigned c2_id,
      unsigned g_idx,
      unsigned edgel_id
      )
{

  dbdet_edgel_chain * c1 = G.fac[c1_id].chain;
  dbdet_edgel_chain * c2 = G.fac[c2_id].chain;

  if (c1->edgels.back()->id == edgel_id)
  {
    if (c2->edgels.front()->id == edgel_id)
      c1->edgels.insert(c1->edgels.end(), ++(c2->edgels.begin()), c2->edgels.end());
    else if(c2->edgels.back()->id == edgel_id)
      c1->edgels.insert(c1->edgels.end(), ++(c2->edgels.rbegin()), c2->edgels.rend());
  }
  else if (c1->edgels.front()->id == edgel_id)
  {
    if (c2->edgels.front()->id == edgel_id)
      c1->edgels.insert(c1->edgels.begin(), c2->edgels.rbegin(), --(c2->edgels.rend()));
    else if(c2->edgels.back()->id == edgel_id)
      c1->edgels.insert(c1->edgels.begin(), c2->edgels.begin(), --(c2->edgels.end()));
  }

  c2->edgels.clear();
  G.var[g_idx].n_facs.clear();
  //G.var[g_idx].merged = true;

  for (vcl_vector<unsigned>::iterator it = G.fac[c1_id].n_vars.begin(); it != G.fac[c1_id].n_vars.end(); it++)
  {
    if(*it == g_idx)
    {
      G.fac[c1_id].n_vars.erase(it);
      break;
    }
  }

  for (vcl_vector<unsigned>::iterator it = G.fac[c2_id].n_vars.begin(); it != G.fac[c2_id].n_vars.end(); it++)
  {
    if(*it == g_idx)
    {
      G.fac[c2_id].n_vars.erase(it);
      break;
    }
  }

  for (unsigned i = 0; i < G.fac[c2_id].n_vars.size(); ++i)
  {
    for (unsigned j = 0; j < G.var[G.fac[c2_id].n_vars[i]].n_facs.size(); ++j)
    {
      if(G.var[G.fac[c2_id].n_vars[i]].n_facs[j] == c2_id)
      {
        G.var[G.fac[c2_id].n_vars[i]].n_facs[j] = c1_id;
        break;
      }
    }
  }

  G.fac[c1_id].n_vars.insert(G.fac[c1_id].n_vars.end(), G.fac[c2_id].n_vars.begin(), G.fac[c2_id].n_vars.end());
  //G.fac[c2_id].removed = true;
}

void dbdet_graphical_model_contour_merge::
compute_texture_hist(
      dbdet_edgel_chain & chain,
      vcl_vector< vnl_vector_fixed<double, 2> > & n,
      y_hist_vector & left, 
      y_hist_vector & right
      )
{
  unsigned const tex_nbr_dist = 3;
  unsigned npts =  chain.edgels.size();


  for (unsigned i = 0; i < y_hist_size; ++i)
    left[i] = right[i] = 0;

  for (unsigned k = 0; k < npts; ++k)
  {
    vgl_point_2d<double> & cur_pt = chain.edgels[k]->pt; 
    for (int l = 1; l <= tex_nbr_dist; ++l)
    {
      unsigned il = vcl_max(0, vcl_min(static_cast<int>(ni()) - 1, static_cast<int>(cur_pt.x() - n[k][0] * l + 0.5)));
      unsigned ir = vcl_max(0, vcl_min(static_cast<int>(ni()) - 1, static_cast<int>(cur_pt.x() + n[k][0] * l + 0.5)));

      unsigned jl = vcl_max(0, vcl_min(static_cast<int>(nj()) - 1, static_cast<int>(cur_pt.y() - n[k][1] * l + 0.5)));
      unsigned jr = vcl_max(0, vcl_min(static_cast<int>(nj()) - 1, static_cast<int>(cur_pt.y() + n[k][1] * l + 0.5)));

      left[tmap_(il, jl)]++;
      right[tmap_(ir, jr)]++;
    }
  }

  for (unsigned i = 0; i < y_hist_size; ++i)
  {
    left[i] /= npts * tex_nbr_dist;
    right[i] /= npts * tex_nbr_dist;
  }
}

void dbdet_graphical_model_contour_merge::
deep_copy_cfg(
      dbdet_curve_fragment_graph & CFG,
      dbdet_curve_fragment_graph & newCFG
      )
{
  newCFG.clear();
  newCFG.resize(CFG.cFrags.size());
  for (dbdet_edgel_chain_list_const_iter it=CFG.frags.begin(); it != CFG.frags.end(); it++)
    newCFG.insert_fragment(new dbdet_edgel_chain(*(*it)));
}
