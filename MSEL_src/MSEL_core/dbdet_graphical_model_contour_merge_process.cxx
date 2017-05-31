// This is brcv/seg/dbdet/pro/dbdet_graphical_model_contour_merge_process.cxx

//:
// \file

#include "dbdet_graphical_model_contour_merge_process.h"


#include "vidpro1_image_storage.h"
#include "vidpro1_vsol2D_storage.h"
#include "vidpro1_vsol2D_storage_sptr.h"

#include "dbdet_edgemap_storage.h"
#include "dbdet_edgemap_storage_sptr.h"
#include "dbdet_sel_storage.h"
#include "dbdet_sel_storage_sptr.h"

#include <vcl_vector.h>
#include <vcl_string.h>
#include <vil/vil_image_resource.h>
#include <vil/vil_new.h>
#include <vil/vil_image_view.h>

#include <vsol/vsol_point_2d.h>
#include <vsol/vsol_point_2d_sptr.h>
#include <vsol/vsol_line_2d.h>
#include <vsol/vsol_line_2d_sptr.h>

#include "dbdet_edgemap_sptr.h"
//#include "dbdet_third_order_edge_det.h"
#include <vil/vil_fill.h>
#include <vil/vil_border.h>

#include "dbtest_root_dir.h"
#include <vnl/vnl_matrix.h>
#include <vil/vil_rgb.h>
#include <vil/vil_convert.h>
#include "dbdet_yuliang_features.h"
#include "dbdet_graphical_model_contour_merge.h"


//: Constructor
dbdet_graphical_model_contour_merge_process::dbdet_graphical_model_contour_merge_process()
{

  vcl_string root = dbtest_root_dir();
  vcl_string base_path = root + "/MSEL_core/test_data";


  fb = new dbdet_filter_bank(base_path);
  tex = new dbdet_texton_classifier((base_path + "/tex/tex.txt").c_str());
/*
   0.0000000e+00   1.2100564e-01   1.1583408e-01   8.5717724e-02   1.5083448e-01   1.1946507e-01   2.2051502e+00   7.6632618e-01   9.1834344e-01
   1.0000000e+00   1.1751265e-01   1.0735821e-01   1.3995593e-01   1.4768553e-01   1.0967336e-01   1.6652428e+00   3.4109466e-01   4.1954811e-01
  -1.7082151e-01  -2.4023619e-01  -2.8814772e-01  -2.8576244e-01   4.9989011e-02  -1.5149038e-01  -4.1761749e-01   9.2396348e-01   1.5311652e-02
*/

  if( !parameters()->add( "fmean[0]"   , "-fmean_0" , 0.0000000e+00) ||
      !parameters()->add( "fmean[1]"   , "-fmean_1" , 1.2100564e-01) ||
      !parameters()->add( "fmean[2]"   , "-fmean_2" , 1.1583408e-01) ||
      !parameters()->add( "fmean[3]"   , "-fmean_3" , 8.5717724e-02) ||
      !parameters()->add( "fmean[4]"   , "-fmean_4" , 1.5083448e-01) ||
      !parameters()->add( "fmean[5]"   , "-fmean_5" , 1.1946507e-01) ||
      !parameters()->add( "fmean[6]"   , "-fmean_6" , 2.2051502e+00) ||
      !parameters()->add( "fmean[7]"   , "-fmean_7" , 7.6632618e-01) ||
      !parameters()->add( "fmean[8]"   , "-fmean_8" , 9.1834344e-01) ||
      !parameters()->add( "fstd[0]"   , "-fstd_0" , 1.0000000e+00) ||
      !parameters()->add( "fstd[1]"   , "-fstd_1" , 1.1751265e-01) ||
      !parameters()->add( "fstd[2]"   , "-fstd_2" , 1.0735821e-01) ||
      !parameters()->add( "fstd[3]"   , "-fstd_3" , 1.3995593e-01) ||
      !parameters()->add( "fstd[4]"   , "-fstd_4" , 1.4768553e-01) ||
      !parameters()->add( "fstd[5]"   , "-fstd_5" , 1.0967336e-01) ||
      !parameters()->add( "fstd[6]"   , "-fstd_6" , 1.6652428e+00) ||
      !parameters()->add( "fstd[7]"   , "-fstd_7" , 3.4109466e-01) ||
      !parameters()->add( "fstd[8]"   , "-fstd_8" , 4.1954811e-01) ||
      !parameters()->add( "beta[0]"   , "-beta_0" , -1.7658682e-01) ||
      !parameters()->add( "beta[1]"   , "-beta_1" , -2.4023619e-01) ||
      !parameters()->add( "beta[2]"   , "-beta_2" , -2.8814772e-01) ||
      !parameters()->add( "beta[3]"   , "-beta_3" , -2.8576244e-01) ||
      !parameters()->add( "beta[4]"   , "-beta_4" , 4.9989011e-02) ||
      !parameters()->add( "beta[5]"   , "-beta_5" , -1.5149038e-01) ||
      !parameters()->add( "beta[6]"   , "-beta_6" , -4.1761749e-01) ||
      !parameters()->add( "beta[7]"   , "-beta_7" , 9.2396348e-01) ||
      !parameters()->add( "beta[8]"   , "-beta_8" , 1.5311652e-02) ||
      !parameters()->add( "fmean_g[0]"   , "-fmean_g_0" , 0.0000000e+00) ||
      !parameters()->add( "fmean_g[1]"   , "-fmean_g_1" , 7.6632618e-01) ||
      !parameters()->add( "fstd_g[0]"   , "-fstd_g_0" , 1.0000000e+00) ||
      !parameters()->add( "fstd_g[1]"   , "-fstd_g_1" , 3.4109466e-01) ||
      !parameters()->add( "beta_g[0]"   , "-beta_g_0" , -1.7658682e-01) ||
      !parameters()->add( "beta_g[1]"   , "-beta_g_1" , 1.0618483e+00)
    )
  {
    vcl_cerr << "ERROR: Adding parameters in " __FILE__ << vcl_endl;
  }
}

//: Destructor
dbdet_graphical_model_contour_merge_process::~dbdet_graphical_model_contour_merge_process()
{
  delete fb;
  delete tex;
}


//: Clone the process
bpro1_process*
dbdet_graphical_model_contour_merge_process::clone() const
{
  return new dbdet_graphical_model_contour_merge_process(*this);
}


//: Return the name of this process
vcl_string
dbdet_graphical_model_contour_merge_process::name()
{
  return "Graphical Model Contour Merge";
}


//: Return the number of input frame for this process
int
dbdet_graphical_model_contour_merge_process::input_frames()
{
  return 3;
}


//: Return the number of output frames for this process
int
dbdet_graphical_model_contour_merge_process::output_frames()
{
  return 1;
}

//: Provide a vector of required input types
vcl_vector< vcl_string > dbdet_graphical_model_contour_merge_process::get_input_type()
{
  vcl_vector< vcl_string > to_return;
  to_return.push_back( "image" );
  to_return.push_back( "edge_map" );
  to_return.push_back( "sel" );
  return to_return;
}


//: Provide a vector of output types
vcl_vector< vcl_string > dbdet_graphical_model_contour_merge_process::get_output_type()
{
  vcl_vector<vcl_string > to_return;
  to_return.push_back( "sel" );
  return to_return;
}


//: Execute the process
bool
dbdet_graphical_model_contour_merge_process::execute()
{
  if ( input_data_.size() != 3 ){
    vcl_cout << "In dbdet_graphical_model_contour_merge_process::execute() - not exactly three"
             << " inputs\n";
    return false;
  }
  clear_output();

  vcl_cout << "Graphical model contour merge...\n";
  vcl_cout.flush();

  // get image from the storage class
  vidpro1_image_storage_sptr frame_image;
  frame_image.vertical_cast(input_data_[0][0]);
  vil_image_resource_sptr image_sptr = frame_image->get_image();
  vil_image_view<vil_rgb<vxl_byte> > image_view = image_sptr->get_view(0, image_sptr->ni(), 0, image_sptr->nj() );

  dbdet_sel_storage_sptr input_sel;
  input_sel.vertical_cast(input_data_[0][2]);
  dbdet_curve_fragment_graph& CFG = input_sel->CFG();

  vcl_cout << "Input #fragments: " << CFG.frags.size() << vcl_endl;

  dbdet_edgemap_storage_sptr input_edgemap;
  input_edgemap.vertical_cast(input_data_[0][1]);
  dbdet_edgemap_sptr EM = input_edgemap->get_edgemap();

  dbdet_sel_storage_sptr output_sel = dbdet_sel_storage_new();
  output_sel->set_EM(input_sel->EM());
  dbdet_curve_fragment_graph &newCFG = output_sel->CFG();

  //get the parameters
  y_params_0_vector fmean_sem, fstd_sem, beta_sem;
  y_params_1_vector fmean_geom, fstd_geom, beta_geom;

  parameters()->get_value( "-fmean_0", fmean_sem[0]);
  parameters()->get_value( "-fmean_1", fmean_sem[1]);
  parameters()->get_value( "-fmean_2", fmean_sem[2]);
  parameters()->get_value( "-fmean_3", fmean_sem[3]);
  parameters()->get_value( "-fmean_4", fmean_sem[4]);
  parameters()->get_value( "-fmean_5", fmean_sem[5]);
  parameters()->get_value( "-fmean_6", fmean_sem[6]);
  parameters()->get_value( "-fmean_7", fmean_sem[7]);
  parameters()->get_value( "-fmean_8", fmean_sem[8]);

  parameters()->get_value( "-fstd_0", fstd_sem[0]);
  parameters()->get_value( "-fstd_1", fstd_sem[1]);
  parameters()->get_value( "-fstd_2", fstd_sem[2]);
  parameters()->get_value( "-fstd_3", fstd_sem[3]);
  parameters()->get_value( "-fstd_4", fstd_sem[4]);
  parameters()->get_value( "-fstd_5", fstd_sem[5]);
  parameters()->get_value( "-fstd_6", fstd_sem[6]);
  parameters()->get_value( "-fstd_7", fstd_sem[7]);
  parameters()->get_value( "-fstd_8", fstd_sem[8]);

  parameters()->get_value( "-beta_0", beta_sem[0]);
  parameters()->get_value( "-beta_1", beta_sem[1]);
  parameters()->get_value( "-beta_2", beta_sem[2]);
  parameters()->get_value( "-beta_3", beta_sem[3]);
  parameters()->get_value( "-beta_4", beta_sem[4]);
  parameters()->get_value( "-beta_5", beta_sem[5]);
  parameters()->get_value( "-beta_6", beta_sem[6]);
  parameters()->get_value( "-beta_7", beta_sem[7]);
  parameters()->get_value( "-beta_8", beta_sem[8]);

  parameters()->get_value( "-fmean_g_0", fmean_geom[0]);
  parameters()->get_value( "-fmean_g_1", fmean_geom[1]);
  parameters()->get_value( "-fstd_g_0", fstd_geom[0]);
  parameters()->get_value( "-fstd_g_1", fstd_geom[1]);
  parameters()->get_value( "-beta_g_0", beta_geom[0]);
  parameters()->get_value( "-beta_g_1", beta_geom[1]);

  for (unsigned i = 0; i < y_params_0_size; ++i)
    beta_sem[i] /= fstd_sem[i];

  for (unsigned i = 0; i < y_params_1_size; ++i)
    beta_geom[i] /= fstd_geom[i];

  vcl_vector<vil_image_view<double> > decomposed = fb->decompose(image_view);
	vnl_matrix<unsigned> tmap = tex->classify(decomposed);

  vcl_cout << "Merging contours...\n"; 
  // perfrom third-order edge detection with these parameters
  dbdet_graphical_model_contour_merge cm(image_view, *EM, tmap);

  
  cm.dbdet_merge_contour(CFG, beta_geom, fmean_geom, beta_sem, fmean_sem, newCFG);

  // create the output storage class
  
  output_data_[0].push_back(output_sel);

  vcl_cout << "Output #fragments: " << newCFG.frags.size() << vcl_endl;
  vcl_cout << "done!" << vcl_endl;
  vcl_cout.flush();

  return true;
}

bool
dbdet_graphical_model_contour_merge_process::finish()
{
  return true;
}
