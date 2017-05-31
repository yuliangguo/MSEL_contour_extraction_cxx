// This is brcv/seg/dbdet/pro/dbdet_load_edg_process_h_

#ifndef dbdet_load_edg_process_h_
#define dbdet_load_edg_process_h_

//:
// \file
// \brief A process for loading a .EDG file into the current frame
// \author Amir Tamrakar
// \date 06/06/04
//
//
// \verbatim
//  Modifications
// \endverbatim

#include "bpro1_process.h"
#include "bpro1_parameters.h"
#include "vidpro1_vsol2D_storage.h"
#include "vidpro1_vsol2D_storage_sptr.h"
#include <vcl_vector.h>
#include "dbdet_edgemap_sptr.h"

//: This process loads a .EDG file into the current frame

class dbdet_load_edg_process : public bpro1_process
{
public:
  dbdet_load_edg_process();
  ~dbdet_load_edg_process() {}

  //: Clone the process
  virtual bpro1_process* clone() const;
  
  vcl_string name();
  vcl_vector< vcl_string > get_input_type();
  vcl_vector< vcl_string > get_output_type();
  
  int input_frames() {
    return 1;
  }
  int output_frames() {
    return num_frames_;
  }
  
  bool execute();
  bool finish() {
    return true;
  }
  
  bool loadEDG(vcl_string filename);


protected:
  int num_frames_;


};

#endif // dbdet_load_edg_process_h_
