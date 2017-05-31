// This is brcv/seg/dbdet/pro/dbdet_save_cem_process.h

#ifndef dbdet_save_cem_process_h_
#define dbdet_save_cem_process_h_

//:
// \file
// \brief A process for saving contours(edgel chains) as .cem files
// \author Amir Tamrakar
// \date 11/15/07
//
// \verbatim
//  Modifications
// \endverbatim

#include "bpro1_process.h"
#include "bpro1_parameters.h"

#include <vcl_vector.h>

//: This process saves a contour fragment graph into a .cem file
class dbdet_save_cem_process : public bpro1_process
{
public:
  dbdet_save_cem_process();
  ~dbdet_save_cem_process() {}

  //: Clone the process
  virtual bpro1_process* clone() const;
  
  vcl_string name();
  vcl_vector< vcl_string > get_input_type();
  vcl_vector< vcl_string > get_output_type();
  
  int input_frames() {
    return 1;
  }
  int output_frames() {
    return 1;
  }
  
  bool execute();
  bool finish() {
    return true;
  }
  
protected:
};

#endif // dbdet_save_cem_process_h_
