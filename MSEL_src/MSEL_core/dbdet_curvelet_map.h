// This is dbdet_curvelet_map.h
#ifndef dbdet_curvelet_map_h
#define dbdet_curvelet_map_h
//:
//\file
//\brief Curvelet Map data structure
//\author Amir Tamrakar
//\date 11/25/07
//
//\verbatim
//  Modifications
//    Amir Tamrakar  Instead of storing the curvelets with the edgels, I've created this 
//                   data structure to hold them so that the edgel classes are lighter and
//                   more relevant to other algorithms.
//
//    Amir Tamrakar  Moved the parameters of the curvelet formation into a separate params class
//
//\endverbatim

#include <vcl_vector.h>
#include <vcl_list.h>

#include "dbdet_edgel.h"
#include "dbdet_curvelet.h"
#include "dbdet_edgemap_sptr.h"

//: curvelet formation parameters
class dbdet_curvelet_params
{
public:
  dbdet_curve_model::curve_type C_type; ///< curve model type

  double rad_;    ///< radius of the grouping neighborhood around each edgel
  double dtheta_; ///< assumed uncertainty in orientation
  double dpos_;   ///< assumed uncertainty in position
  bool badap_uncer_; ///< get the uncertainty from the edges
  double token_len_; ///< Length of the edgel token (puts bounds on legal curvature)
  double max_k_;  ///< maximum curvature of curves in the curve bundle
  double max_gamma_; ///< maximum curvature derivative of curves in the curve bundle

  unsigned nrad_; ///< the range of edgel cells to search (this should contain all edgels within rad_) 
  double gap_;
  unsigned maxN_; ///< largest curvelet size to form
  bool centered_; ///< anchored centered curvelet
  bool bidirectional_; ///< form curvelets in both directions from the anchored

  dbdet_curvelet_params(dbdet_curve_model::curve_type ctype=dbdet_curve_model::CC2, double rad=7.0, double gap=1.0,double dtheta=15, 
                        double dpos=0.5, bool adap_uncer=true, double token_len=0.7, double max_k=0.5, double max_gamma=0.05,
                        bool centered=true, bool bidir=false): 
    C_type(ctype), rad_(rad),gap_(gap),dtheta_(dtheta), dpos_(dpos), badap_uncer_(adap_uncer), token_len_(token_len), 
    max_k_(max_k), max_gamma_(max_gamma), nrad_((unsigned) vcl_ceil(rad)+1), maxN_(2*nrad_), 
    centered_(centered), bidirectional_(bidir)
  {}

  ~dbdet_curvelet_params(){}
};

//: This class stores the map of curvelets formed by the SEL edge linker.
//  This is an intermediate structure before actual edge linking occurs.
class dbdet_curvelet_map
{
public:
  //: The edgemap on which these curvelets have been formed 
  //  (due to this smart pointer to the edgemap, the curvelets remain valid even if the edgemap is deleted elsewhere)
  dbdet_edgemap_sptr EM_; 

  //: various parameters used for forming this map
  dbdet_curvelet_params params_;

  //: The curvelet map, indexed by edgel IDs
  vcl_vector<cvlet_list > map_;

  //: The curvelet map for the other direction (only for DHT mode)
  vcl_vector<cvlet_list > map2_;

  //: constructor
  dbdet_curvelet_map(dbdet_edgemap_sptr EM=0, dbdet_curvelet_params params=dbdet_curvelet_params());

  //: destructor
  ~dbdet_curvelet_map();

  //: to check if the Curvelet map is valid
  bool is_valid() { return map_.size()>0 && map_.size()==EM_->num_edgels(); }

  //: set the edgemap
  void set_edgemap(dbdet_edgemap_sptr EM) { EM_ = EM; resize(EM->num_edgels()); }

  //: set the linking parameters
  void set_parameters(dbdet_curvelet_params params) { params_ = params; }

  //: access the curvelets for an edge using its id
  const cvlet_list& curvelets(unsigned id) const { return map_[id]; }
  cvlet_list& curvelets(unsigned id) { return map_[id]; }

  cvlet_list& Rcurvelets(unsigned id) { return map2_[id]; }

  //: resize the graph
  void resize(unsigned size);

  //: clear the graph
  void clear();

  //: clear all the curvelets in the graph
  void clear_all_curvelets();

  //: add a curvelet to an edgel
  void add_curvelet(dbdet_curvelet* curvelet, bool dir=true);

  //: remove a curvelet from this edgel
  void remove_curvelet(dbdet_curvelet* curvelet);

  //: delete all the curvelets formed by this edgel
  void delete_all_curvelets(dbdet_edgel* e);

  //: does this curvelet exist at this edgel?
  dbdet_curvelet* does_curvelet_exist(dbdet_edgel* e, vcl_deque<dbdet_edgel*> & chain);

  //: does the given pair exist on the ref edgel?
  dbdet_curvelet* find_pair(dbdet_edgel* ref, dbdet_edgel* eA, dbdet_edgel* eB);

  //: does the given triplet exist on the ref edgel?
  dbdet_curvelet* find_triplet(dbdet_edgel* ref, dbdet_edgel* eA, dbdet_edgel* eB, dbdet_edgel* eC);

  //: return largest curvelet formed by an edgel
  dbdet_curvelet* largest_curvelet(dbdet_edgel* e);

  friend class dbdet_edge_map;
};

#endif // dbdet_curvelet_map_h
