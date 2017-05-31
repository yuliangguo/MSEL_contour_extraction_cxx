// This is brcv/seg/dbdet/pro/dbdet_graphical_model_contour_merge_process.h
#ifndef dbdet_graphical_model_contour_merge_process_h_
#define dbdet_graphical_model_contour_merge_process_h_

#include "bpro1_process.h"
#include "bpro1_parameters.h"

#include "dbdet_filter_util.h"
#include "dbdet_filter_bank.h"
#include "dbdet_texton_classifier.h"

//: Third order edge detector
class dbdet_graphical_model_contour_merge_process : public bpro1_process 
{
  dbdet_filter_bank * fb;
  dbdet_texton_classifier * tex;
public:

  dbdet_graphical_model_contour_merge_process();
  virtual ~dbdet_graphical_model_contour_merge_process();

  //: Clone the process
  virtual bpro1_process* clone() const;

  vcl_string name();

  int input_frames();
  int output_frames();


  vcl_vector< vcl_string > get_input_type();
  vcl_vector< vcl_string > get_output_type();

  bool execute();
  bool finish();

protected:
  int num_frames_;

};

#endif
