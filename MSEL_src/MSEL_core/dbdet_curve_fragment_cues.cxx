#include "dbdet_curve_fragment_cues.h"
#include "bgld_diffgeom.h"
#include <vil/algo/vil_colour_space.h>
#include <vcl_cmath.h>
#include <vcl_iostream.h>

void dbdet_curve_fragment_cues::
compute_all_cues(
      const dbdet_edgel_chain &c, 
      y_feature_vector *features_ptr // indexed by the enum
      )
{
  unsigned const npts = c.edgels.size();
  const dbdet_edgel_list e = c.edgels;
  //const vil_image_view<vxl_uint_32> dt = *dt_;
  y_feature_vector &features = *features_ptr;
  features[y_features::Y_ONE] = 1;

  cuvature_cues(c, features_ptr);
  hsv_gradient_cues(c, features_ptr);
  features[y_features::Y_EDGE_SPARSITY] = lateral_edge_sparsity_cue(c, features_ptr);
  //mean_conf = mean(cfrag(:,4));

  // compute average edge strength (mean_conf)

  double conf=0;
  for (dbdet_edgel_list_const_iter eit=c.edgels.begin(); eit != c.edgels.end(); eit++) {
    assert(*eit);
    conf += (*eit)->strength;
  }
  if (npts)
    features[y_features::Y_MEAN_CONF] = conf / npts;
}

void
dbdet_curve_fragment_cues::
cuvature_cues(
      const dbdet_edgel_chain &c, 
      y_feature_vector *features_ptr // indexed by the enum
      )
{
  unsigned const npts = c.edgels.size();
  if (npts < 2) // curvature is only computed for two samples
    return;
  y_feature_vector &features = *features_ptr;
  // curvature
  vnl_vector<double> k;
  vcl_vector<vgl_point_2d<double> > points;
  points.reserve(npts);

  //to vector of points..
  for (unsigned i = 0; i < npts; ++i)
    points.push_back(c.edgels[i]->pt);

  bgld_compute_curvature(points, &k);
  assert(k.size() == npts);

  for (unsigned i=0; i < npts; ++i)
    if (vnl_math::isnan(k[i]))
      k[i] = 0;

  { // wiggliness is the time curvature change sign
    features[y_features::Y_WIGG] = 0;
    for (unsigned i=0; i + 1 < npts; ++i)
    {
      double val = k[i+1] * k[i];
      if (val < 0.0 && vnl_math::abs(val) > epsilon)
      {
        features[y_features::Y_WIGG] += i + 1;
      }
    }
    features[y_features::Y_WIGG] /= npts;
  }
 
  features[y_features::Y_ABS_K] = k.one_norm() / npts;
}

void
dbdet_curve_fragment_cues::
hsv_gradient_cues(
      const dbdet_edgel_chain &c, 
      y_feature_vector *features_ptr // indexed by the enum
    )
{
  // examine local_dist neighborhood around curve along normals
  unsigned const npts = c.edgels.size();
  if (!npts) return;
  vcl_vector< vnl_vector_fixed<double, 2> > n;
  n.reserve(npts);
  vcl_vector<vgl_point_2d<double> > points;
  points.reserve(npts);

  //to vector of points..
  for (unsigned i = 0; i < npts; ++i)
    points.push_back(c.edgels[i]->pt);
  bgld_compute_normals(points, &n);

  // get neighborhood points to be examined
  y_feature_vector &features = *features_ptr;

  features[y_features::Y_HUE_GRAD] = 0;
  features[y_features::Y_SAT_GRAD] = 0;
  features[y_features::Y_BG_GRAD]  = 0;

  for (unsigned i=0; i < npts; ++i) {
    unsigned left_x  = static_cast<unsigned>(points[i].x() - local_dist_ * n[i][0] + 0.5);
    unsigned left_y  = static_cast<unsigned>(points[i].y() - local_dist_ * n[i][1] + 0.5);
    unsigned right_x = static_cast<unsigned>(points[i].x() + local_dist_ * n[i][0] + 0.5);
    unsigned right_y = static_cast<unsigned>(points[i].y() + local_dist_ * n[i][1] + 0.5);

    double hue_left, hue_right,
           sat_left, sat_right,
           bg_left, bg_right;

    vil_rgb<vxl_byte> rgb = im(left_x,left_y);
    vil_colour_space_RGB_to_HSV<double>(rgb.r, rgb.g, rgb.b, 
        &hue_left, &sat_left, &bg_left);
   rgb = im(right_x,right_y);
   vil_colour_space_RGB_to_HSV<double>(rgb.r, rgb.g, rgb.b, 
        &hue_right, &sat_right, &bg_right);

    // TODO test if imae indexing is (x,y) or (y,x)
    features[y_features::Y_SAT_GRAD] += /*vcl_abs*/(sat_left - sat_right);
    features[y_features::Y_BG_GRAD]  += /*vcl_abs*/(bg_left - bg_right) / 255.;
    // TODO need to make angle difference to make sense
    features[y_features::Y_HUE_GRAD] += /*vcl_abs*/(hue_left - hue_right)/360.;
  }
  features[y_features::Y_HUE_GRAD] /= npts;
  features[y_features::Y_SAT_GRAD] /= npts;
  features[y_features::Y_BG_GRAD]  /= npts;
}

double dbdet_curve_fragment_cues::
lateral_edge_sparsity_cue(
    const dbdet_edgel_chain &c,
    y_feature_vector *features_ptr // indexed by the enum
    )
{
  y_feature_vector &features = *features_ptr;
  unsigned const npts = c.edgels.size();
  unsigned total_edges = 0;
  /*if (use_dt()) {
    for (unsigned i=0; i < npts; ++i) {
      // locate bucket
      unsigned p_i = static_cast<unsigned>(e[i]->pt.x()+0.5);
      unsigned p_j = static_cast<unsigned>(e[i]->pt.y()+0.5);
      
      if (dt(p_i, p_j) > nbr_width_)
        continue;

      // visit all nearby edgels and count the number within a distance
      // mark already visited edgels

      // TODO:optimize access to be row-first
      for (int d_i = -nbr_width_; di < nbr_width_; ++d_i) {
        for (int d_j = -nbr_width_; dj < nbr_width_; ++d_j) {
          if (not_visited(p_i + d_i, p_j + d_j)) {
            unsigned nh_x = static_cast<unsigned>(p_i + d_i);
            unsigned nh_y = static_cast<unsigned>(p_j + d_j);
            total_edges += em.cell(nh_x,nh_y).size();
            mark_visited(nh_x,nh_y);
          }
        }
      }
    }
  } else { // no dt
  }
  */
  assert(!use_dt());

  const dbdet_edgel_list &e = c.edgels;
  /*for (unsigned i=0; i < npts; ++i) {
    unsigned px = static_cast<unsigned>(e[i]->pt.x()+0.5);
    unsigned py = static_cast<unsigned>(e[i]->pt.y()+0.5);

    assert (px < ni());
    assert (py < nj());
    mark_visited(px, py);
  }

  for (unsigned i=0; i < npts; ++i) {
    // locate bucket
    unsigned p_i = static_cast<unsigned>(e[i]->pt.x()+0.5);
    unsigned p_j = static_cast<unsigned>(e[i]->pt.y()+0.5);
    
    // visit all nearby edgels and count the number within a distance
    // mark already visited edgels

    // TODO:optimize access to be row-first
    for (int d_i = -nbr_width_; d_i < nbr_width_; ++d_i) {
      for (int d_j = -nbr_width_; d_j < nbr_width_; ++d_j) {
        if (not_visited(p_i + d_i, p_j + d_j)) {
          unsigned nh_x = static_cast<unsigned>(p_i + d_i);
          unsigned nh_y = static_cast<unsigned>(p_j + d_j);
          total_edges += (em_.cell(nh_x,nh_y).size() > 0 ? 1 : 0);
          mark_visited(nh_x,nh_y);
        }
      }
    }
  }
  visited_id_++;*/
  

  //naive implementation
  int w = ni();
  int h = nj();

  for(int i = 0; i < w; ++i)
    for(int j = 0; j < h; ++j)
      mask(i, j) = 0;

  for(int k = 0; k < npts; ++k)
  {
    int px = static_cast<int>(e[k]->pt.x()+0.5);
    int py = static_cast<int>(e[k]->pt.y()+0.5);
    for(int i = vcl_max(px - static_cast<int>(nbr_width_), 0); i < vcl_min(w, 1 + px + static_cast<int>(nbr_width_)); i++)
      for(int j = vcl_max(py - static_cast<int>(nbr_width_), 0); j < vcl_min(h, 1 + py + static_cast<int>(nbr_width_)); j++)
        mask(i, j) = 1;
  }

  for(int i = 0; i < npts; ++i)
  {
    unsigned px = static_cast<unsigned>(e[i]->pt.x()+0.5);
    for(int j = 0; j < npts; ++j)
    {
      unsigned py = static_cast<unsigned>(e[j]->pt.y()+0.5);
      mask(px, py) = 0;
    }
  }

  for(int i = 0; i < w; ++i)
  {
    for(int j = 0; j < h; ++j)
    {
      (em_.cell(i, j).size() > 0 && mask(i, j)) ? total_edges++ : 0;
    }
  }

  features[y_features::Y_LEN] = euclidean_length(c);
  return total_edges / (features[y_features::Y_LEN] == 0 ? 1.0 : features[y_features::Y_LEN]);
}
