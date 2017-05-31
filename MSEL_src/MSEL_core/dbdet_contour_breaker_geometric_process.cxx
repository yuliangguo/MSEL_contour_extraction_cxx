// This is brcv/seg/dbdet/pro/dbdet_contour_breaker_geometric_process.cxx

//:
// \file

#include "dbdet_contour_breaker_geometric_process.h"


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

//#include <dbtest/dbtest_root_dir.h>
#include <vnl/vnl_matrix.h>
#include <vil/vil_rgb.h>
#include <vil/vil_convert.h>
#include "dbdet_yuliang_features.h"
#include "dbdet_contour_breaker.h"

//: Constructor
dbdet_contour_breaker_geometric_process::dbdet_contour_breaker_geometric_process()
{
  if( !parameters()->add( "fmean[0]"   , "-fmean_0" , 0.0000000e+00) ||
      !parameters()->add( "fmean[1]"   , "-fmean_1" , 7.6632618e-01) ||
      !parameters()->add( "fstd[0]"   , "-fstd_0" , 1.0000000e+00) ||
      !parameters()->add( "fstd[1]"   , "-fstd_1" , 3.4109466e-01) ||
      !parameters()->add( "beta[0]"   , "-beta_0" , -1.7658682e-01) ||
      !parameters()->add( "beta[1]"   , "-beta_1" , 1.0618483e+00)
    )
  {
    vcl_cerr << "ERROR: Adding parameters in " __FILE__ << vcl_endl;
  }
}

//: Destructor
dbdet_contour_breaker_geometric_process::~dbdet_contour_breaker_geometric_process()
{
}


//: Clone the process
bpro1_process*
dbdet_contour_breaker_geometric_process::clone() const
{
  return new dbdet_contour_breaker_geometric_process(*this);
}


//: Return the name of this process
vcl_string
dbdet_contour_breaker_geometric_process::name()
{
  return "Geometric Contour Breaker";
}


//: Return the number of input frame for this process
int
dbdet_contour_breaker_geometric_process::input_frames()
{
  return 3;
}


//: Return the number of output frames for this process
int
dbdet_contour_breaker_geometric_process::output_frames()
{
  return 1;
}

//: Provide a vector of required input types
vcl_vector< vcl_string > dbdet_contour_breaker_geometric_process::get_input_type()
{
  vcl_vector< vcl_string > to_return;
  to_return.push_back( "image" );
  to_return.push_back( "edge_map" );
  to_return.push_back( "sel" );
  return to_return;
}


//: Provide a vector of output types
vcl_vector< vcl_string > dbdet_contour_breaker_geometric_process::get_output_type()
{
  vcl_vector<vcl_string > to_return;
  to_return.push_back( "sel" );
  return to_return;
}


//: Execute the process
bool
dbdet_contour_breaker_geometric_process::execute()
{
  if ( input_data_.size() != 3 ){
    vcl_cout << "In dbdet_contour_breaker_geometric_process::execute() - not exactly three"
             << " inputs\n";
    return false;
  }
  clear_output(1);

  vcl_cout << "Gemometric contour breaker...\n";
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
  y_params_1_vector fmean_geom, fstd_geom, beta_geom;

  parameters()->get_value( "-fmean_0", fmean_geom[0]);
  parameters()->get_value( "-fmean_1", fmean_geom[1]);
  parameters()->get_value( "-fstd_0", fstd_geom[0]);
  parameters()->get_value( "-fstd_1", fstd_geom[1]);
  parameters()->get_value( "-beta_0", beta_geom[0]);
  parameters()->get_value( "-beta_1", beta_geom[1]);

  for (unsigned i = 0; i < y_params_1_size; ++i)
    beta_geom[i] /= fstd_geom[i];

  //vcl_vector<vil_image_view<double> > decomposed = fb->decompose(image_view);
	//vnl_matrix<unsigned> tmap = tex->classify(decomposed);
  vnl_matrix<unsigned> tmap(image_view.ni(), image_view.nj());
  // perfrom third-order edge detection with these parameters
  dbdet_contour_breaker cb(image_view, *EM, tmap);

  
  cb.dbdet_contour_breaker_geom(CFG, beta_geom, fmean_geom, newCFG);

  // create the output storage class
  
  output_data_[0].push_back(output_sel);

  vcl_cout << "Output #fragments: " << newCFG.frags.size() << vcl_endl;
  vcl_cout << "done!" << vcl_endl;
  vcl_cout.flush();

  return true;
}

bool
dbdet_contour_breaker_geometric_process::finish()
{
  return true;
}
