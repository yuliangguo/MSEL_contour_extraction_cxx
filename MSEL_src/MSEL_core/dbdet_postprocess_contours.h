// This is brcv/seg/dbdet/algo/dbdet_postprocess_contours.h
#ifndef dbdet_postprocess_contours_h
#define dbdet_postprocess_contours_h
//:
//\file
//\brief Some methods to postprocess contours
//\author Amir Tamrakar
//\date 11/25/07
//
//\verbatim
//  Modifications
//   Amir Tamrakar Moved these functions out of the edge linking algorithms 
//                 Because they are common to multiple algorithms
//\endverbatim

#include "dbdet_curve_fragment_graph.h"

//********************************************************************//
// Functions for post processing of Linked Contours
//********************************************************************//

//Functions for computing various contour properties
double avg_strength(dbdet_edgel_chain* chain);
double mean_contrast(dbdet_edgel_chain* chain);
double left_app_std(dbdet_edgel_chain* chain);
double right_app_std(dbdet_edgel_chain* chain);

double avg_d2f(dbdet_edgel_chain* chain);

void compute_curvatures(dbdet_edgel_chain* chain, vcl_vector<double> &ks);
double avg_curvature(dbdet_edgel_chain* chain);
double max_curvature(dbdet_edgel_chain* chain);

//: Prune contours on the basis of the length of the curve
void prune_contours_by_length(dbdet_curve_fragment_graph& curve_frag_graph, 
                              double len_thresh=10.0);

//: Prune contours on the basis of the length of the curve
void prune_contours_by_strength(dbdet_curve_fragment_graph& curve_frag_graph, 
                              double strength_thresh=3.0);

//: Prune contours on the basis of average contrast and length of the curve
void prune_contours(dbdet_curve_fragment_graph& curve_frag_graph, 
                    double len_thresh=0.0, 
                    double strength_thresh=0.0, 
                    double contrast_thresh=0.0, 
                    double adap_thresh_fac=0.0, 
                    double d2f_thresh=0.0, 
                    double k_thresh=0.0);

//: just testing: use a moving window
void appearance_based_post_processing(dbdet_curve_fragment_graph& curve_frag_graph, unsigned win_len, double adap_thresh);

//: post process to break contours at high curvature points
void post_process_based_on_curvature(dbdet_curve_fragment_graph& curve_frag_graph, double k_thresh=1.0);

//: Post process to link distinct but connected curve fragments
void post_process_to_link_fragments(dbdet_curve_fragment_graph& curve_frag_graph);

//: Perform post-processing to remove ambiuous segments between the ends of contour fragments
void post_process_to_remove_ambiguity(dbdet_curve_fragment_graph& curve_frag_graph);


#endif // dbdet_postprocess_contours_h
