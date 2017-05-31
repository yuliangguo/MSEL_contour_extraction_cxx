// This is brcv/seg/dbdet/pro/dbdet_contour_ranker_process.h
#ifndef dbdet_contour_ranker_process_h_
#define dbdet_contour_ranker_process_h_

#include "bpro1_process.h"
#include "bpro1_parameters.h"

//: Third order edge detector
class dbdet_contour_ranker_process : public bpro1_process 
{
public:

  dbdet_contour_ranker_process();
  virtual ~dbdet_contour_ranker_process();

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
