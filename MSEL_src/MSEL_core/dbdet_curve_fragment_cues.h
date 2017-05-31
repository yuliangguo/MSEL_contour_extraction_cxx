#include <vnl/vnl_vector_fixed.h>
#include <vbl/vbl_array_2d.h>
#include <vxl_config.h>
#include <vcl_limits.h>
#include <vil/vil_image_view.h>
#include <vil/vil_rgb.h>
#include <vil/vil_border.h>
#include "dbdet_edgel.h"
#include "dbdet_edgemap.h"
#include "dbdet_yuliang_features.h"


// This is dbdet_curve_fragment_cues.h
#ifndef dbdet_curve_fragment_cues_h
#define dbdet_curve_fragment_cues_h
//:
//\file
//\brief Compute several geometric and appearance cues of a curve fragment
//\author Ricardo Fabbri (rfabbri), Brown University  (rfabbri.github.io)
//\date 05/25/2016 21:57:40 BRT
//

//static const vxl_uint_32 dbdet_curve_fragment_cues_unvisited = vcl_numeric_limits<vxl_uint_32>::max();

//: Compute curve fragment cues for many curves.
// Holds state information such as distance transform and auxiliary buffers,
// so that computation os fast for many curve fragments in the same image.
class dbdet_curve_fragment_cues {
public:

  // Make sure the input parameters stay valid
  // while this class is in use.
  dbdet_curve_fragment_cues(
    const vil_image_view<vil_rgb<vxl_byte> > &img,
    const dbdet_edgemap &em
    )
    :
    //visited_img_(img.ni(), img.nj()),
    //visited_id_(0),
    //visited_(vil_border_create_accessor(visited_img_, vil_border_create_constant(visited_img_, dbdet_curve_fragment_cues_unvisited))),
    img_(img),
    im(vil_border_create_accessor(img_,vil_border_create_geodesic(img_))),
    mask(img.ni(), img.nj()),
    //    dt_(dt),
    em_(em),
    dt_(NULL)
  {
    /*if (img.nplanes() != 3){
      std::cerr << "Input must be 3-plane RGB images!" << std::endl;
      abort();
    }*/

    //visited_img_.fill(dbdet_curve_fragment_cues_unvisited);
    // outside indices return false (visited)
    //visited_ = vil_border_create_accessor(visited_img_, vil_border_create_constant(visited_img_, dbdet_curve_fragment_cues_unvisited));
  }

  // Try to speedup lateral edge sparsity computation for many fragments,
  // using with distance transform.
  //
  // Pass in the curve fragment list which will be queried.
  //  void use_dt(const dbdet_edgel_chain_list &frags)
  //  {
    // remove edgels of curves from curve maps
    // perform DT of the remaining map
    // const vil_image_view<vxl_uint_32> &dt
  //    dt_ = &dt;
  //  }

  void compute_all_cues(
      const dbdet_edgel_chain &c, 
      y_feature_vector *features_ptr // indexed by the enum
      );
  void cuvature_cues(
        const dbdet_edgel_chain &c, 
        y_feature_vector *features_ptr // indexed by the enum
        );
  void hsv_gradient_cues(
        const dbdet_edgel_chain &c, 
        y_feature_vector *features_ptr // indexed by the enum
      );
  double lateral_edge_sparsity_cue(
      const dbdet_edgel_chain &c,
      y_feature_vector *features_ptr // indexed by the enum
    );

  static double euclidean_length(const dbdet_edgel_chain &c) {
    double len=0;
    for (unsigned i=0; i+1 < c.edgels.size(); ++i)
      len += vgl_distance(c.edgels[i]->pt, c.edgels[i+1]->pt);
    return len;
  }
private:
  unsigned ni() const { return em_.ncols(); }
  unsigned nj() const { return em_.nrows(); }
  bool use_dt() const { return dt_ != NULL; }
  //bool visited(int i, int j) const { return visited_(i,j) == visited_id_; }
  //bool not_visited(int i, int j) const { return visited_(i,j) != visited_id_; }
  //void mark_visited(int i, int j) { visited_img_(i,j) = visited_id_; }
  const vil_image_view<vil_rgb<vxl_byte> > &img_; // color RGB image
  // make sure image clamps to within bounds
  vil_border_accessor<vil_image_view<vil_rgb<vxl_byte> > > im;
  vil_image_view<vxl_uint_32> *dt_;
  // visited(i,j) = c marks pixels(i,j) as visited by curve c's nhood tube
  // visited(i,j) = UIHNT_MAX marks pixels(i,j) as not visited
  //vil_image_view<vxl_uint_32> visited_img_;
  //vil_border_accessor<vil_image_view<vxl_uint_32> > visited_;
  // each compute_cues() run increments the id to mark traversals in scrap buffer visited_
  // that way we reuse the buffer without clearing it
  //vxl_uint_32 visited_id_;
  vbl_array_2d<bool> mask;
  const dbdet_edgemap &em_;
  static unsigned const local_dist_ = 2; // distance used for local sampling
  static unsigned const nbr_width_ = 3;  // distance used for lateral edge sparsity
  static double const epsilon = 1e-10;
};


#endif // dbdet_curve_fragment_cues_h
