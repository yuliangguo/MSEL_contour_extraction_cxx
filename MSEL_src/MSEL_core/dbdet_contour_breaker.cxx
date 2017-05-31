#include "dbdet_contour_breaker.h"
#include "dbdet_edgel.h"
#include <vcl_algorithm.h>
#include <vcl_limits.h>
#include "bgld_diffgeom.h"
#include <vil/algo/vil_colour_space.h>
#include <vil/vil_border.h>
#include <vcl_iostream.h>
#include <vcl_cassert.h>

void dbdet_contour_breaker::
dbdet_contour_breaker_geom(
      dbdet_curve_fragment_graph & CFG,
      y_params_1_vector & beta,
      y_params_1_vector & fmean,
      dbdet_curve_fragment_graph & newCFG
      )
{
  int const ref_tabel_nbr_range = 2;
  unsigned const not_assigned = -1;

  deep_copy_cfg(CFG, newCFG);

  vcl_vector <dbdet_edgel_chain*> frags(newCFG.frags.size());
  {
    unsigned i = 0;
    for (dbdet_edgel_chain_list_const_iter it=newCFG.frags.begin(); it != newCFG.frags.end(); it++, ++i)
      frags[i] = (*it);
  }

  vcl_vector<double> clen(frags.size(), 0.0);

  double min_contour_len = dbdet_yuliang_const::dbdet_yuliang_const::nbr_len_th * diag_ratio;

  for (unsigned i = 0; i < ni(); ++i)
  {
    for(unsigned j = 0; j < nj(); ++j)
    {
      ref_start_pts(i, j) = not_assigned;
      ref_end_pts(i, j) = not_assigned;
    }
  }

  for (unsigned i = 0; i < frags.size(); ++i)
  {
    clen[i] = euclidean_length(*frags[i]);

    dbdet_edgel & start = *(frags[i]->edgels.front());
    dbdet_edgel & end = *(frags[i]->edgels.back());

    if (start.pt == end.pt && start.deriv == end.deriv)
      continue;

    if (clen[i] > min_contour_len && frags[i]->edgels.size() > dbdet_yuliang_const::nbr_len_th)
    {
      unsigned xi, xf, yi, yf;
      int x, y;
      x = static_cast<int>(start.pt.x() + 0.5);
      y = static_cast<int>(start.pt.y() + 0.5);
      xi = static_cast<unsigned>(vcl_max(x - ref_tabel_nbr_range, 0));
      xf = static_cast<unsigned>(vcl_min(x + ref_tabel_nbr_range + 1, static_cast<int>(ni())));
      yi = static_cast<unsigned>(vcl_max(y - ref_tabel_nbr_range, 0));
      yf = static_cast<unsigned>(vcl_min(y + ref_tabel_nbr_range + 1, static_cast<int>(nj())));
      
      for (unsigned l = xi; l < xf; ++l)
        for (unsigned m = yi; m < yf; ++m)
          ref_start_pts(l, m) = i;

      x = static_cast<int>(end.pt.x() + 0.5);
      y = static_cast<int>(end.pt.y() + 0.5);
      xi = static_cast<unsigned>(vcl_max(x - ref_tabel_nbr_range, 0));
      xf = static_cast<unsigned>(vcl_min(x + ref_tabel_nbr_range + 1, static_cast<int>(ni())));
      yi = static_cast<unsigned>(vcl_max(y - ref_tabel_nbr_range, 0));
      yf = static_cast<unsigned>(vcl_min(y + ref_tabel_nbr_range + 1, static_cast<int>(nj())));
      
      for (unsigned l = xi; l < xf; ++l)
        for (unsigned m = yi; m < yf; ++m)
          ref_end_pts(l, m) = i;
    }  
  }

  for (unsigned i = 0; i < frags.size(); ++i)
  {
    dbdet_edgel & start = *(frags[i]->edgels.front());
    dbdet_edgel & end = *(frags[i]->edgels.back());

    if (start.pt == end.pt && start.deriv == end.deriv)
      continue;

    unsigned npts = frags[i]->edgels.size();

    if(clen[i] > (10 * diag_ratio) && npts > nbr_num_edges + 1)
    {
      vcl_vector<unsigned> cur_break_e_id;
      vcl_vector <int> start_ids(npts);
      vcl_set <unsigned> unique_start_ids;
      vcl_vector <int> end_ids(npts);
      vcl_set <unsigned> unique_end_ids;

      dbdet_edgel_list edgels = frags[i]->edgels;

      for (unsigned k = 0; k < npts; ++k)
      {
        int xx = static_cast<int>(edgels[k]->pt.x() + 0.5);
        int yy = static_cast<int>(edgels[k]->pt.y() + 0.5);
        unsigned x = vcl_max(vcl_min(xx, static_cast<int>(ni()) - 1), 0);
        unsigned y = vcl_max(vcl_min(yy, static_cast<int>(nj()) - 1), 0);
          
        start_ids[k] = ref_start_pts[x][y];
        end_ids[k] = ref_end_pts[x][y];

        if (ref_start_pts[x][y] != not_assigned && ref_start_pts[x][y] != i)
          unique_start_ids.insert(ref_start_pts[x][y]);

        if (ref_end_pts[x][y] != not_assigned && ref_end_pts[x][y] != i)
          unique_end_ids.insert(ref_end_pts[x][y]);
      }

      compute_break_point(frags, i, start_ids, unique_start_ids, true, cur_break_e_id);
      compute_break_point(frags, i, end_ids, unique_end_ids, false, cur_break_e_id);

      if (cur_break_e_id.size() > 0)
      {
        vcl_sort(cur_break_e_id.begin(), cur_break_e_id.end());
        dbdet_edgel_chain copy = *frags[i];
        dbdet_edgel_chain & ref = *frags[i];

        dbdet_edgel_chain * newChain;

        ref.edgels = dbdet_edgel_list(copy.edgels.begin() + cur_break_e_id[cur_break_e_id.size() - 1], copy.edgels.end());
        clen[i] = euclidean_length(ref);

        newChain = new dbdet_edgel_chain();
        newChain->edgels = dbdet_edgel_list(copy.edgels.begin(), copy.edgels.begin() + cur_break_e_id[0] + 1);    
        newCFG.frags.push_back(newChain);
        clen.push_back(euclidean_length(*newChain));
        frags.push_back(newChain);

        unsigned xi, xf, yi, yf;
        int x, y;
        x = static_cast<int>(start.pt.x() + 0.5);
        y = static_cast<int>(start.pt.y() + 0.5);
        xi = static_cast<unsigned>(vcl_max(x - ref_tabel_nbr_range, 0));
        xf = static_cast<unsigned>(vcl_min(x + ref_tabel_nbr_range, static_cast<int>(ni())));
        yi = static_cast<unsigned>(vcl_max(y - ref_tabel_nbr_range, 0));
        yf = static_cast<unsigned>(vcl_min(y + ref_tabel_nbr_range, static_cast<int>(nj())));
        
        for (unsigned l = xi; l < xf; ++l)
          for (unsigned m = yi; m < yf; ++m)
            ref_start_pts(l, m) = newCFG.frags.size()-1;

        for (unsigned k = 1; k < cur_break_e_id.size(); ++k)
        {
          newChain = new dbdet_edgel_chain();
          newChain->edgels = dbdet_edgel_list(copy.edgels.begin() + cur_break_e_id[k-1], copy.edgels.begin() + cur_break_e_id[k] + 1);
          newCFG.frags.push_back(newChain);
          clen.push_back(euclidean_length(*newChain));
          frags.push_back(newChain);
        }   
      }
    }  
  }

  for(unsigned i = 0; i < max_it; ++i)
  {
    unsigned end_j = newCFG.frags.size();
    unsigned j = 0;
 
    for(dbdet_edgel_chain_list_iter it = newCFG.frags.begin(); (it != newCFG.frags.end() && j < end_j); it++,++j)
    {
      dbdet_edgel & start = *((*it)->edgels.front());
      dbdet_edgel & end = *((*it)->edgels.back());
      
      if (start.pt == end.pt && start.deriv == end.deriv)
        continue;

      if (clen[j] > (10 * diag_ratio) && (*it)->edgels.size() > nbr_num_edges / (i + 1.0))
      {
        vcl_vector<double> prob;
        compute_merge_probability_geom(*(*it), vcl_ceil(nbr_num_edges /(2.0 * (i + 1.0))) , beta, fmean, prob);
        double min = vcl_numeric_limits<double>::max();
        double lastVal = min;
        unsigned firstId = -1, lastId = -1;
        for (int k = 0; k < prob.size(); ++k)
        {
          if (prob[k] < min)
          {
            min = prob[k];
            firstId = lastId = k; 
          }
          else if(prob[k] == min)
          {
            lastId = k;
          }
        }

        if (min < dbdet_yuliang_const::merge_th_geom)
        {
          unsigned id = (firstId + lastId + 1) / 2;
          dbdet_edgel_chain copy = *(*it);
          dbdet_edgel_chain & ref = *(*it);

          dbdet_edgel_chain * newChain;

          ref.edgels = dbdet_edgel_list(copy.edgels.begin() + id, copy.edgels.end());
          clen[i] = euclidean_length(ref);

          newChain = new dbdet_edgel_chain();
          newChain->edgels = dbdet_edgel_list(copy.edgels.begin(), copy.edgels.begin() + id + 1);
          newCFG.frags.push_back(newChain);
          clen.push_back(euclidean_length(*newChain));
        }
      }
    }
  }
}

void dbdet_contour_breaker::
compute_break_point(
      vcl_vector<dbdet_edgel_chain*> & frags,
      unsigned frag_id,
      vcl_vector<int> & ids,
      vcl_set<unsigned> & unique_ids,
      bool front,
      vcl_vector<unsigned> & break_e_ids)
{
      unsigned prev_id = 0;
      for (vcl_set<unsigned>::iterator set_it = unique_ids.begin(); set_it != unique_ids.end(); ++set_it)
      {
        dbdet_edgel_chain & chain = (*frags[frag_id]);
        double min = vcl_numeric_limits<double>::max();
        unsigned e_id = -1;

        vgl_point_2d<double> set_pt;
        if(front)
        {
          set_pt = (*frags[(*set_it)]).edgels.front()->pt;
        }
        else
        {
          set_pt = (*frags[(*set_it)]).edgels.back()->pt;
        }

        for (int k = 0; k < ids.size(); ++k)
        {
          if((*set_it) == ids[k])
          {
            double dist = vgl_distance(chain.edgels[k]->pt, set_pt);
            if (dist < min)
            {
              min = dist;
              e_id = k;
            }
          }
        }

        if (e_id < dbdet_yuliang_const::dbdet_yuliang_const::nbr_len_th - 1 || e_id > chain.edgels.size() - dbdet_yuliang_const::dbdet_yuliang_const::nbr_len_th || (e_id - prev_id) < dbdet_yuliang_const::dbdet_yuliang_const::nbr_len_th)
          continue;

        if(min == 0)
        {
          //junction_pts = [junction_pts; cur_c(e_id, :)];
          //introduced_num_junction_points = introduced_num_junction_points + 1; 
          break_e_ids.push_back(e_id);
          continue;
        }

        vgl_vector_2d<double> c_ori = chain.edgels[e_id + 1]->pt - chain.edgels[e_id - 1]->pt;

        unsigned a_id1, a_id2;
        if (front)
        {
          a_id1 = vcl_min(dbdet_yuliang_const::dbdet_yuliang_const::nbr_len_th, static_cast<unsigned>(frags[(*set_it)]->edgels.size() - 1));
          a_id2 = 0;
        }
        else
        {
          a_id1 = frags[(*set_it)]->edgels.size() - 1;
          a_id2 = static_cast<unsigned>(vcl_max(static_cast<int>(frags[(*set_it)]->edgels.size()) - 1 - static_cast<int>(dbdet_yuliang_const::nbr_len_th), 0));
        }
        vgl_vector_2d<double> a_ori = (*frags[(*set_it)]).edgels[a_id1]->pt - (*frags[(*set_it)]).edgels[a_id2]->pt;

        double cos_diff = (c_ori.x() * a_ori.x() + c_ori.y() * a_ori.y()) / 
        (vcl_sqrt(c_ori.x() * c_ori.x() + c_ori.y() * c_ori.y()) * vcl_sqrt(a_ori.x() * a_ori.x() + a_ori.y() * a_ori.y()));

        if(cos_diff > vcl_cos(vnl_math::pi / 6.0))
          continue;

        break_e_ids.push_back(e_id);
        prev_id = e_id;
      }
}

void dbdet_contour_breaker::
compute_merge_probability_geom(
      dbdet_edgel_chain & chain,
      unsigned nbr_range_th,
      y_params_1_vector & beta,
      y_params_1_vector & fmean,
      vcl_vector<double> & prob
      )
{
  nbr_range_th = vcl_min(nbr_range_th, dbdet_yuliang_const::nbr_len_th);
  prob.resize(chain.edgels.size());
  
  for (unsigned i = 0; i < prob.size(); ++i)
    prob[i] = 1.0;

  for (unsigned i = nbr_range_th; i < prob.size() - nbr_range_th; ++i)
  {
    vgl_vector_2d<double> a = chain.edgels[i]->pt - chain.edgels[i-nbr_range_th]->pt;
    vgl_vector_2d<double> b = chain.edgels[i+nbr_range_th]->pt - chain.edgels[i]->pt;

    double cos_diff = (a.x() * b.x() + a.y() * b.y()) / 
                      (vcl_sqrt(a.x() * a.x() + a.y() * a.y()) * vcl_sqrt(b.x() * b.x() + b.y() * b.y()));
    //prob(i) = 1 / (1 + exp(-([1 geom_diff] - fmean_1)*beta_1'));
    prob[i] =  1.0 / (1.0 + exp(-((1.0 - fmean[0]) * beta[0] + (cos_diff - fmean[1]) * beta[1])));
  }
}

void dbdet_contour_breaker::
dbdet_contour_breaker_semantic(
      dbdet_curve_fragment_graph & CFG,
      y_params_0_vector & beta,
      y_params_0_vector & fmean,
      dbdet_curve_fragment_graph & newCFG
      )
{
  //dbdet_curve_fragment_graph newCFG(CFG);
  deep_copy_cfg(CFG, newCFG);

  vcl_vector<double> clen(newCFG.frags.size(), 0.0);

  for(unsigned i = 0; i < max_it; ++i)
  {
    unsigned end_j = newCFG.frags.size();
    unsigned j = 0;

    for(dbdet_edgel_chain_list_iter it = newCFG.frags.begin(); (it != newCFG.frags.end() && j < end_j); it++, ++j)
    {
      dbdet_edgel & start = *((*it)->edgels.front());
      dbdet_edgel & end = *((*it)->edgels.back());


      clen[j] = (clen[j] != 0.0 ? clen[j] : euclidean_length(*(*it)));

      if (start.pt == end.pt && start.deriv == end.deriv)
        continue;


      if (clen[j] > (30 * diag_ratio) && (*it)->edgels.size() > (3 * nbr_num_edges) / i)
      {
        vcl_vector<double> prob;
        compute_merge_probability_semantic(*(*it), static_cast<unsigned>(nbr_num_edges /(2 * i) + 0.5) , beta, fmean, prob);
        double min = vcl_numeric_limits<double>::max();
        double lastVal = min;
        unsigned firstId = -1, lastId = -1;
        for (int k = 0; k < prob.size(); ++k)
        {
          if (prob[k] < min)
          {
            min = prob[k];
            firstId = lastId = k; 
          }
          else if(prob[k] == min)
          {
            lastId = k;
          }
        }
        if (min < dbdet_yuliang_const::merge_th_sem)
        {
          unsigned id = (firstId + lastId + 1) / 2;
          dbdet_edgel_chain copy = *(*it);
          dbdet_edgel_chain & ref = *(*it);

          dbdet_edgel_chain * newChain;

          ref.edgels = dbdet_edgel_list(copy.edgels.begin() + id, copy.edgels.end());
          clen[j] = euclidean_length(ref);

          newChain = new dbdet_edgel_chain();
          newChain->edgels = dbdet_edgel_list(copy.edgels.begin(), copy.edgels.begin() + id + 1);
          newCFG.frags.push_back(newChain);
          clen.push_back(euclidean_length(*newChain));
        }
      }
    }
  }
}

void dbdet_contour_breaker::
compute_merge_probability_semantic(
      dbdet_edgel_chain & chain,/*hsv_img, edge_map, tmap,*/
      unsigned nbr_range_th,
      y_params_0_vector & beta,
      y_params_0_vector & fmean,
      vcl_vector<double> & prob
      )
{
  unsigned const local_dist = 1;
  unsigned const nbr_width = 3;
  
  unsigned npts =  chain.edgels.size();

  vnl_vector<double> k;
  vcl_vector< vnl_vector_fixed<double, 2> > n;
  vcl_vector<vgl_point_2d<double> > points;
  points.reserve(npts);

  prob.resize(npts);
  for (unsigned i = 0; i < npts; ++i)
    prob[i] = 1.0;

  //to vector of points..
  for (unsigned i = 0; i < npts; ++i)
    points.push_back(chain.edgels[i]->pt);

  bgld_compute_curvature(points, &k);
  bgld_compute_normals(points, &n);

  for (unsigned i=0; i < npts; ++i)
    if (vnl_math::isnan(k[i]))
      k[i] = 0;

  vcl_vector<double> k_cum(npts, 0.0), hue_cum(npts, 0.0), sat_cum(npts, 0.0), bg_cum(npts, 0.0), wigg_cum(npts, 0.0);
  double last_k, last_hue, last_sat, last_bg, last_wigg;
  last_k = last_hue = last_sat = last_bg = last_wigg = 0;

  // make sure image clamps to within bounds
  vil_border_accessor<vil_image_view<vil_rgb<vxl_byte> > >
      im = vil_border_create_accessor(img_,vil_border_create_geodesic(img_));

  for (unsigned i=0; i < npts; ++i) {
    unsigned left_x  = static_cast<unsigned>(points[i].x() - local_dist * n[i][0] + 0.5);
    unsigned left_y  = static_cast<unsigned>(points[i].y() - local_dist * n[i][1] + 0.5);
    unsigned right_x = static_cast<unsigned>(points[i].x() + local_dist * n[i][0] + 0.5);
    unsigned right_y = static_cast<unsigned>(points[i].y() + local_dist * n[i][1] + 0.5);

    

    double hue_left, hue_right,
           sat_left, sat_right,
           bg_left, bg_right;

    vil_rgb<vxl_byte> rgb = im(left_x,left_y);
    vil_colour_space_RGB_to_HSV<double>(rgb.r, rgb.g, rgb.b, 
        &hue_left, &sat_left, &bg_left);
    rgb = im(right_x,right_y);
    vil_colour_space_RGB_to_HSV<double>(rgb.r, rgb.g, rgb.b, 
        &hue_right, &sat_right, &bg_right);

    
    last_sat += (sat_left - sat_right);
    last_bg += (bg_left - bg_right) / 255.;
    last_hue += (hue_left - hue_right)/360.;
    
    hue_cum[i] = last_hue;
    sat_cum[i] = last_sat;
    bg_cum[i] = last_bg;

    last_k += vcl_abs(k[i]);
    k_cum[i] = last_k;
    
    double sign = k[i] * (i+1 < npts ? k[i+1] : 0.0);
    last_wigg += (sign < 0 && vcl_abs(sign) > dbdet_yuliang_const::epsilon) ? 1.0 : 0.0;
    wigg_cum[i] = last_wigg;
  }

  vcl_vector<double> edge_sparsity_cum(npts, 0.0);
  compute_edge_sparsity_integral(chain, n, nbr_width, edge_sparsity_cum);
/**
[~, ~, texton_hist_left_cum, texton_hist_right_cum] = compute_texture_hist_integral(x,y,N, nbr_width, tmap);
**/

  vcl_vector<y_hist_vector> texton_hist_left, texton_hist_right;
  compute_texture_hist_integral(chain, n, nbr_width, texton_hist_left, texton_hist_right);
  y_params_0_vector features;
  features[y_params_0::Y_ONE] = 1.0;

  for (unsigned i = nbr_range_th; i < npts - nbr_range_th; ++i)
  {
    features[y_params_0::Y_BG_GRAD] = (2.0 * bg_cum[i] - bg_cum[i-nbr_range_th] - bg_cum[i+nbr_range_th]) / nbr_range_th;
    features[y_params_0::Y_SAT_GRAD] = (2.0 * sat_cum[i] - sat_cum[i-nbr_range_th] - sat_cum[i+nbr_range_th]) / nbr_range_th;
    features[y_params_0::Y_HUE_GRAD] = (2.0 * hue_cum[i] - hue_cum[i-nbr_range_th] - hue_cum[i+nbr_range_th]) / nbr_range_th;
    features[y_params_0::Y_ABS_K] = (2.0 * k_cum[i] - k_cum[i-nbr_range_th] - k_cum[i+nbr_range_th]) / nbr_range_th;
    features[y_params_0::Y_EDGE_SPARSITY] = (2.0 * edge_sparsity_cum[i] - edge_sparsity_cum[i-nbr_range_th] - edge_sparsity_cum[i+nbr_range_th]) / nbr_range_th;
    features[y_params_0::Y_WIGG] = (2.0 * wigg_cum[i] - wigg_cum[i-nbr_range_th] - wigg_cum[i+nbr_range_th]) / nbr_range_th;

   /**
    hist_1_left = (texton_hist_left_cum(:,i) - texton_hist_left_cum(:,i-nbr_range_th))/nbr_range_th;
    hist_2_left = (texton_hist_left_cum(:,i+nbr_range_th) - texton_hist_left_cum(:, i))/nbr_range_th;
    left_dist = pdist2(hist_1_left',hist_2_left','chisq');
    
    hist_1_right = (texton_hist_right_cum(:,i) - texton_hist_right_cum(:,i-nbr_range_th))/nbr_range_th;
    hist_2_right = (texton_hist_right_cum(:,i+nbr_range_th) - texton_hist_right_cum(:, i))/nbr_range_th;
    right_dist = pdist2(hist_1_right',hist_2_right','chisq');    
    
    texture_diff = left_dist + right_dist;

     'chisq'
%   The chi-squared distance between two vectors is defined as:
%    d(x,y) = sum( (xi-yi)^2 / (xi+yi) ) / 2;
%   The chi-squared distance is useful when comparing histograms.
    **/
    
    //texture_diff
    features[y_params_0::Y_TEXTURE] = 0.0;
    for (unsigned k = 0; k < y_hist_size; ++k)
    {
      double v1 = texton_hist_left[i][k] - texton_hist_left[i-nbr_range_th][k];
      double v2 = texton_hist_left[i][k+nbr_range_th] - texton_hist_left[i][k];
      features[y_params_0::Y_TEXTURE] += (v1 - v2) * (v1 - v2) / (v1 + v2);

      v1 = texton_hist_right[i][k] - texton_hist_right[i-nbr_range_th][k];
      v2 = texton_hist_right[i][k+nbr_range_th] - texton_hist_right[i][k];
      features[y_params_0::Y_TEXTURE] += (v1 - v2) * (v1 - v2) / (v1 + v2);
    }
    features[y_params_0::Y_TEXTURE] /= (2 * nbr_range_th);

    vgl_vector_2d<double> a = chain.edgels[i]->pt - chain.edgels[i-nbr_range_th]->pt;
    vgl_vector_2d<double> b = chain.edgels[i+nbr_range_th]->pt - chain.edgels[i]->pt;

    //geom_diff
    features[y_params_0::Y_GEOM] = (a.x() * b.x() + a.y() * b.y()) / 
                      (vcl_sqrt(a.x() * a.x() + a.y() * a.y()) * vcl_sqrt(b.x() * b.x() + b.y() * b.y()));
    double sum = 0.0;
    for (unsigned j = 0; j < y_params_0_size; ++j)
    {
      sum += (features[j] - fmean[j]) * beta[j];
    }
    prob[i] = 1.0 / (1.0 + exp(-sum));
  }
}

void dbdet_contour_breaker::
compute_edge_sparsity_integral(
      dbdet_edgel_chain & chain,
      vcl_vector< vnl_vector_fixed<double, 2> > & n,
      unsigned nbr_width,
      vcl_vector<double> & edge_sparsity
      )
{
  unsigned npts =  chain.edgels.size();
  edge_sparsity.resize(npts);
  double last = 0.0;

  for (unsigned k = 0; k < npts; ++k)
  {
    vgl_point_2d<double> & cur_pt = chain.edgels[k]->pt; 
    unsigned xi = static_cast<unsigned>(vcl_max(0, vcl_min(static_cast<int>(ni()), static_cast<int>(cur_pt.x() - n[k][0] * nbr_width + 0.5))));
    unsigned xf = static_cast<unsigned>(vcl_max(0, vcl_min(static_cast<int>(ni()), static_cast<int>(cur_pt.x() + n[k][0] * nbr_width + 0.5))));
    unsigned yi = static_cast<unsigned>(vcl_max(0, vcl_min(static_cast<int>(nj()), static_cast<int>(cur_pt.x() - n[k][1] * nbr_width + 0.5))));
    unsigned yf = static_cast<unsigned>(vcl_max(0, vcl_min(static_cast<int>(nj()), static_cast<int>(cur_pt.x() + n[k][1] * nbr_width + 0.5))));
    if(xi > xf) vcl_swap(xi, xf);
    if(yi > yf) vcl_swap(yi, yf);

    double sum = 0;
    for (unsigned i = xi; i < xf; ++i)
    {
      for (unsigned j = yi; j < yf; ++j)
      {
        sum = (em_.cell(i, j).size() > 0 ? sum + 1.0: sum);
      }
    }
    last += sum;
    edge_sparsity[k] = last;
  }
}

void dbdet_contour_breaker::
compute_texture_hist_integral(
      dbdet_edgel_chain & chain,
      vcl_vector< vnl_vector_fixed<double, 2> > & n,
      unsigned nbr_width,
      vcl_vector<y_hist_vector> & texton_hist_left,
      vcl_vector<y_hist_vector> & texton_hist_right
      )
{
  unsigned npts =  chain.edgels.size();
  texton_hist_left.resize(npts);
  texton_hist_right.resize(npts);
  double nbr_width_inv = 1.0 / nbr_width;
  
  y_hist_vector last_left = texton_hist_left[0], last_right = texton_hist_right[0];

  for (unsigned i = 0; i < y_hist_size; ++i)
    last_left[i] = last_right[i] = 0;

  for (unsigned k = 0; k < npts; ++k)
  {

    vgl_point_2d<double> & cur_pt = chain.edgels[k]->pt;

    int xl = vcl_max(0, vcl_min(static_cast<int>(ni()) - 1, static_cast<int>(cur_pt.x() - n[k][0] * nbr_width + 0.5)));
    int xr = vcl_max(0, vcl_min(static_cast<int>(ni()) - 1, static_cast<int>(cur_pt.x() + n[k][0] * nbr_width + 0.5)));
    int yl = vcl_max(0, vcl_min(static_cast<int>(nj()) - 1, static_cast<int>(cur_pt.y() - n[k][1] * nbr_width + 0.5)));
    int yr = vcl_max(0, vcl_min(static_cast<int>(nj()) - 1, static_cast<int>(cur_pt.y() + n[k][1] * nbr_width + 0.5)));
    int x = vcl_max(0, vcl_min(static_cast<int>(ni()) - 1, static_cast<int>(cur_pt.x() + 0.5)));
    int y = vcl_max(0, vcl_min(static_cast<int>(nj()) -1, static_cast<int>(cur_pt.y() + 0.5)));

    texton_hist_left[k] = last_left;
    texton_hist_right[k] = last_right;

    for (int i = 0; i <= nbr_width; ++i)
    {
      unsigned il = static_cast<unsigned>(x + i * (xl - x) * nbr_width_inv + 0.5);
      unsigned ir = static_cast<unsigned>(x + i * (xr - x) * nbr_width_inv + 0.5);

      unsigned jl = static_cast<unsigned>(y + i * (yl - y) * nbr_width_inv + 0.5);
      unsigned jr = static_cast<unsigned>(y + i * (yr - y) * nbr_width_inv + 0.5);

      texton_hist_left[k][tmap_(il, jl)]++;
      texton_hist_right[k][tmap_(ir, jr)]++;
    }
    last_left = texton_hist_left[k];
    last_right = texton_hist_right[k];
  }
}

void dbdet_contour_breaker::
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
