// This is brcv/seg/dbdet/pro/dbdet_contour_ranker_process.cxx

//:
// \file

#include "dbdet_contour_ranker_process.h"


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
#include "dbdet_curve_fragment_ranker.h"
#include <vcl_algorithm.h>

bool my_comp(const vcl_pair<double, unsigned> & a, const vcl_pair<double, unsigned> & b)
{
  return a.first > b.first;
}

//: Constructor
dbdet_contour_ranker_process::dbdet_contour_ranker_process()
{
  if( !parameters()->add( "Best #n contours(0=all)"   , "-nfrags" , 0) ||
      !parameters()->add( "Threshold"   , "-thresh" , 0.0) ||
      !parameters()->add( "Minimum contour length(edgels)"   , "-minlen" , 0) ||
      !parameters()->add( "fmean[0]"   , "-fmean_0" , 0.0000000e+00) ||
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
      !parameters()->add( "beta[8]"   , "-beta_8" , 1.5311652e-02)
    )
  {
    vcl_cerr << "ERROR: Adding parameters in " __FILE__ << vcl_endl;
  }
}

//: Destructor
dbdet_contour_ranker_process::~dbdet_contour_ranker_process()
{
}


//: Clone the process
bpro1_process*
dbdet_contour_ranker_process::clone() const
{
  return new dbdet_contour_ranker_process(*this);
}


//: Return the name of this process
vcl_string
dbdet_contour_ranker_process::name()
{
  return "Contour Ranker";
}


//: Return the number of input frame for this process
int
dbdet_contour_ranker_process::input_frames()
{
  return 3;
}


//: Return the number of output frames for this process
int
dbdet_contour_ranker_process::output_frames()
{
  return 1;
}

//: Provide a vector of required input types
vcl_vector< vcl_string > dbdet_contour_ranker_process::get_input_type()
{
  vcl_vector< vcl_string > to_return;
  to_return.push_back( "image" );
  to_return.push_back( "edge_map" );
  to_return.push_back( "sel" );
  return to_return;
}


//: Provide a vector of output types
vcl_vector< vcl_string > dbdet_contour_ranker_process::get_output_type()
{
  vcl_vector<vcl_string > to_return;
  to_return.push_back( "sel" );
  return to_return;
}


//: Execute the process
bool
dbdet_contour_ranker_process::execute()
{
  if ( input_data_.size() != 3 ){
    vcl_cout << "In dbdet_contour_ranker_process::execute() - not exactly three"
             << " inputs\n";
    return false;
  }
  clear_output(1);

  vcl_cout << "Contour ranker...\n";
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
  y_trained_parameters param;
  double thresh;
  int nfrags, minlen;

  parameters()->get_value( "-nfrags", nfrags);
  parameters()->get_value( "-thresh", thresh);
  parameters()->get_value( "-minlen", minlen);

  parameters()->get_value( "-fmean_0", param[0][0]);
  parameters()->get_value( "-fmean_1", param[0][1]);
  parameters()->get_value( "-fmean_2", param[0][2]);
  parameters()->get_value( "-fmean_3", param[0][3]);
  parameters()->get_value( "-fmean_4", param[0][4]);
  parameters()->get_value( "-fmean_5", param[0][5]);
  parameters()->get_value( "-fmean_6", param[0][6]);
  parameters()->get_value( "-fmean_7", param[0][7]);
  parameters()->get_value( "-fmean_8", param[0][8]);

  parameters()->get_value( "-fstd_0", param[1][0]);
  parameters()->get_value( "-fstd_1", param[1][1]);
  parameters()->get_value( "-fstd_2", param[1][2]);
  parameters()->get_value( "-fstd_3", param[1][3]);
  parameters()->get_value( "-fstd_4", param[1][4]);
  parameters()->get_value( "-fstd_5", param[1][5]);
  parameters()->get_value( "-fstd_6", param[1][6]);
  parameters()->get_value( "-fstd_7", param[1][7]);
  parameters()->get_value( "-fstd_8", param[1][8]);

  parameters()->get_value( "-beta_0", param[2][0]);
  parameters()->get_value( "-beta_1", param[2][1]);
  parameters()->get_value( "-beta_2", param[2][2]);
  parameters()->get_value( "-beta_3", param[2][3]);
  parameters()->get_value( "-beta_4", param[2][4]);
  parameters()->get_value( "-beta_5", param[2][5]);
  parameters()->get_value( "-beta_6", param[2][6]);
  parameters()->get_value( "-beta_7", param[2][7]);
  parameters()->get_value( "-beta_8", param[2][8]);

  vnl_vector<double> rank;
  dbdet_curve_fragment_ranker(CFG.frags, EM, image_view, param, &rank);

  vcl_vector<vcl_pair<double, unsigned> > ri(rank.size());
  for(int i=0; i<ri.size(); ++i)
  {
    ri[i].first = rank[i];
    ri[i].second = i;
  }

  vcl_sort(ri.begin(), ri.end(), my_comp);

  nfrags = nfrags == 0 ? rank.size() : nfrags;

  vcl_vector<bool> disp(rank.size(), false);
  for(int i=0; i< vcl_min(nfrags, (int)rank.size()); ++i)
    disp[ri[i].second] = (ri[i].first >= thresh) ? true : false;

  vcl_cout << "Max rank: (" << ri[0].first << ", " << ri[0].second << ") Min rank: (" << ri[ri.size()-1].first << ", " << ri[ri.size()-1].second << ")" << vcl_endl;

  newCFG.clear();
  newCFG.resize(CFG.cFrags.size());
  dbdet_edgel_chain_list_const_iter it=CFG.frags.begin();

  for(int i=0; (i < disp.size() && it != CFG.frags.end()); ++i, it++)
  {
    if(disp[i] && (*it)->edgels.size() > minlen)
      newCFG.insert_fragment(new dbdet_edgel_chain(*(*it)));
  }
  // create the output storage class
  vcl_cout << "Output #fragments: " << newCFG.frags.size() << vcl_endl;
  output_data_[0].push_back(output_sel);
  vcl_cout << "done!" << vcl_endl;
  vcl_cout.flush();

  return true;
}

bool
dbdet_contour_ranker_process::finish()
{
  return true;
}
