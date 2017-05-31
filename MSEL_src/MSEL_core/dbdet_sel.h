// This is brcv/seg/dbdet/algo/dbdet_sel.h
#ifndef dbdet_sel_h
#define dbdet_sel_h
//:
//\file
//\brief Templatized Symbolic edge linking alogorithm(new version)
//\author Amir Tamrakar
//\date 03/15/06
//
//\verbatim
//  Modifications
//  Amir Tamrakar 09/05/06   Removed all the other classes to their own files
//  Amir Tamrakar 09/07/06   Templatized this class to enable it to operate on
//                           different curve models
//  Amir Tamrakar 09/11/06   Moved the meat of the algorithm to a new base class
//                           and left all the functions that depend on the template
//                           in this class.
//  Ricardo Fabbri 09/08/09  Various optimizations to the linker
//
//\endverbatim

#include "dbdet_sel_utils.h"
#include "dbdet_sel_base.h"

//: A templatized subclass that can work with different curve models
template <class curve_model>
class dbdet_sel : public dbdet_sel_base
{
public:

  //: constructor
  dbdet_sel<curve_model>(dbdet_edgemap_sptr edgemap, 
                         dbdet_curvelet_map& cvlet_map, 
                         dbdet_edgel_link_graph& edge_link_graph, 
                         dbdet_curve_fragment_graph& curve_frag_graph,
                         dbdet_curvelet_params cvlet_params=dbdet_curvelet_params()) :
      dbdet_sel_base(edgemap, cvlet_map, edge_link_graph, curve_frag_graph, cvlet_params)
  {
    curvelet_map_.params_.C_type = curve_model().type; //set the right curvelet type
  }

  //: destructor
  virtual ~dbdet_sel<curve_model>(){}

  //: form a curve hypothesis of the appropriate model given a pair of edgels
  inline curve_model* form_a_hypothesis(dbdet_edgel* ref_e, dbdet_edgel* e2, bool &ref_first, 
      bool forward=true, bool centered=true, bool leading=true)
  {
    // First check for consistency in the appearance information
    if (app_usage_==2){
      if (e2->left_app->dist(*(ref_e->left_app))>app_thresh_ ||
          e2->right_app->dist(*(ref_e->right_app))>app_thresh_)
          return 0; //appearance is not consistent
    }

    // Do the dot product test to determine if this pair of edgels can produce a valid hypothesis
    // The two possible hypotheses are:
    //    a) (ref_e->e2)
    //    b) (e2->ref_e)
    //
    
    const double ref_dir = dbdet_vPointPoint(ref_e->pt, e2->pt); //reference dir (ref_e->e2)

    //ref centered grouping
    if (centered) {
      if (dbdet_dot(ref_dir, ref_e->tangent)>0) {
        if (forward){
          ref_first = true;
          return new curve_model(ref_e, e2, ref_e, dpos_, dtheta_, token_len_, max_k_, max_gamma_, badap_uncer_);
        }
        else {
          ref_first = false;
          return new curve_model(e2, ref_e, ref_e, dpos_, dtheta_, token_len_, max_k_, max_gamma_, badap_uncer_);
        }
      }
      else {
        if (forward){
          ref_first = false;
          return new curve_model(e2, ref_e, ref_e, dpos_, dtheta_, token_len_, max_k_, max_gamma_, badap_uncer_);
        }
        else {
          ref_first = true;
          return new curve_model(ref_e, e2, ref_e, dpos_, dtheta_, token_len_, max_k_, max_gamma_, badap_uncer_);
        }
      }
    }
    else { //not centered
      if (dbdet_dot(ref_dir, ref_e->tangent)>0 && forward && leading) {
        ref_first = true;
        return new curve_model(ref_e, e2, ref_e, dpos_, dtheta_, token_len_, max_k_, max_gamma_, badap_uncer_);
      }
      if (dbdet_dot(ref_dir, ref_e->tangent)<0) {
        if (forward && !leading){
          ref_first = false;
          return new curve_model(e2, ref_e, ref_e, dpos_, dtheta_, token_len_, max_k_, max_gamma_, badap_uncer_);
        }
        if (!forward){
          ref_first = true;
          return new curve_model(ref_e, e2, ref_e, dpos_, dtheta_, token_len_, max_k_, max_gamma_, badap_uncer_);
        }
      }
    }
    return 0;
  }

  virtual void form_an_edgel_pair(dbdet_edgel* ref_e, dbdet_edgel* e2);
  virtual void form_an_edgel_triplet(dbdet_curvelet* /*p1*/, dbdet_curvelet* /*p2*/){}
  virtual void form_an_edgel_quad(dbdet_curvelet* /*t1*/, dbdet_curvelet* /*t2*/){}
  virtual void build_curvelets_greedy_for_edge(dbdet_edgel* eA, unsigned max_size_to_group, 
      bool use_flag=false, bool forward=true,  bool centered=true, bool leading=true);

  //: check if two curvelets are consistent (is there an intersection in their curve bundles?)
  bool are_curvelets_consistent(dbdet_curvelet* cvlet1, dbdet_curvelet* cvlet2)
  {
    curve_model *cm1 = static_cast<curve_model *>(cvlet1->curve_model);
    curve_model *cm2 = static_cast<curve_model *>(cvlet2->curve_model);
    curve_model* new_cm = cm1->intersect(cm2);

    bool consistent = new_cm->bundle_is_valid();
    //house cleaning
    delete new_cm; 
    
    return consistent;
  }

  //: form an edgel grouping from an ordered list of edgemap_->edgels
  virtual dbdet_curvelet* form_an_edgel_grouping(dbdet_edgel* ref_e, vcl_deque<dbdet_edgel*> &edgel_chain, 
      bool forward=true,  bool centered=true, bool leading=true);
};


#endif // dbdet_sel_h
