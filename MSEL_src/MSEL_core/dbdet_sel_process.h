// This is brcv/seg/dbdet/pro/dbdet_sel_process.h
#ifndef dbdet_sel_process_h_
#define dbdet_sel_process_h_

//:
// \file
// \brief Process to perform local edge linking
// \author Amir Tamrakar
// modified by Yuliang Jul/29/2010  can choose not to form link graph
// \date 03/12/06
//
// \verbatim
//  Modifications
// \endverbatim

#include "bpro1_process.h"
#include "bpro1_parameters.h"

//: This process traces contours of binary regions in an image to produce subpixel contours
class dbdet_sel_process : public bpro1_process 
{
public:

  dbdet_sel_process();
  virtual ~dbdet_sel_process();

  //: Clone the process
  virtual bpro1_process* clone() const;

  vcl_string name();

  int input_frames();
  int output_frames();

  vcl_vector< vcl_string > get_input_type();
  vcl_vector< vcl_string > get_output_type();

  bool execute();
  bool finish();

  void get_parameters();

protected:
  //various parameters
  double nrad, gap, dx, dt;
  bool badap_uncer;
  unsigned curve_model_type;
  double token_len, max_k, max_gamma;

  unsigned grouping_algo;
  unsigned max_size_to_group;
  unsigned cvlet_type;
  bool bCentered_grouping;
  bool bBidirectional_grouping;

  bool bFormCompleteCvletMap;
  bool bFormLinkGraph;
  bool b_use_all_cvlets;
  unsigned app_usage;
  double app_thresh;
  
  unsigned linkgraph_algo;
  unsigned min_size_to_link;
  
  unsigned linking_algo;
  unsigned num_link_iters;

  bool bGetfinalcontours;

  bool bmergefrags;

};

#endif
