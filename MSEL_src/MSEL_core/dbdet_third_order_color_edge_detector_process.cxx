// This is brcv/seg/dbdet/pro/dbdet_third_order_color_edge_detector_process.cxx

//:
// \file

#include "dbdet_third_order_color_edge_detector_process.h"

#include "vidpro1_image_storage.h"
#include "vidpro1_vsol2D_storage.h"
#include "vidpro1_vsol2D_storage_sptr.h"

#include "dbdet_edgemap_storage.h"
#include "dbdet_edgemap_storage_sptr.h"
#include <vil/vil_save.h>

#include <vcl_vector.h>
#include <vcl_string.h>
#include <vul/vul_timer.h>
#include <vbl/vbl_array_2d.h>
#include <vnl/vnl_math.h>
#include <vnl/vnl_erf.h>
#include <vil/vil_image_resource.h>
#include <vil/vil_new.h>
#include <vil/vil_image_view.h>
#include <vil/vil_convert.h>

#include <brip/brip_vil_float_ops.h>
#include <bil/algo/bil_color_conversions.h>

#include <vsol/vsol_point_2d.h>
#include <vsol/vsol_point_2d_sptr.h>
#include <vsol/vsol_line_2d.h>
#include <vsol/vsol_line_2d_sptr.h>
#include <vil/vil_fill.h>
#include <vil/vil_border.h>

#include "dbdet_third_order_color_detector.h"

//: Constructor
dbdet_third_order_color_edge_detector_process::dbdet_third_order_color_edge_detector_process()
{
  vcl_vector<vcl_string> gradient_operator_choices;
  gradient_operator_choices.push_back("Gaussian");       //0
  gradient_operator_choices.push_back("h0-operator");    //1
  gradient_operator_choices.push_back("h1-operator");    //2
  
  vcl_vector<vcl_string> convolution_choices;
  convolution_choices.push_back("2-D");            //0
  convolution_choices.push_back("1-D");            //1

  vcl_vector<vcl_string> color_conversion_choices;
  color_conversion_choices.push_back("Use RGB");    //0
  color_conversion_choices.push_back("Use IHS");    //1
  color_conversion_choices.push_back("Use Lab");    //2
  color_conversion_choices.push_back("Use Luv");    //3

  vcl_vector<vcl_string> parabola_fit_type;
  parabola_fit_type.push_back("3-point fit");      //0
  parabola_fit_type.push_back("9-point fit");      //1

  if( !parameters()->add( "Load Component Images"    , "-bLoadComps"   , false ) ||
      !parameters()->add( "Color Space conversion"   , "-col_conv" , color_conversion_choices, 2) ||
      !parameters()->add( "Gradient Operator"   , "-grad_op" , gradient_operator_choices, 0) ||
      !parameters()->add( "Convolution Algo:"   , "-conv_algo" , convolution_choices, 0) ||
      !parameters()->add( "Sigma (Gaussian)"    , "-sigma"   , 1.0 ) ||
      !parameters()->add( "Gradient Magnitude Threshold"   , "-thresh" , 2.0 ) ||
      !parameters()->add( "Interpolation factor [2^N], N= "  , "-int_factor" , 1 ) ||
      !parameters()->add( "Parabola Fit type"   , "-parabola_fit" , parabola_fit_type, 0) ||
      !parameters()->add( "Reduce edgel tokens "  , "-breduce" , false ))
  {
    vcl_cerr << "ERROR: Adding parameters in " __FILE__ << vcl_endl;
  }
}


//: Destructor
dbdet_third_order_color_edge_detector_process::~dbdet_third_order_color_edge_detector_process()
{
}


//: Clone the process
bpro1_process*
dbdet_third_order_color_edge_detector_process::clone() const
{
  return new dbdet_third_order_color_edge_detector_process(*this);
}


//: Return the name of this process
vcl_string
dbdet_third_order_color_edge_detector_process::name()
{
  return "Third Order Color Edge Detector";
}


//: Return the number of input frame for this process
int
dbdet_third_order_color_edge_detector_process::input_frames()
{
  return 1;
}


//: Return the number of output frames for this process
int
dbdet_third_order_color_edge_detector_process::output_frames()
{
  return 1;
}


//: Provide a vector of required input types
vcl_vector< vcl_string > dbdet_third_order_color_edge_detector_process::get_input_type()
{
  vcl_vector< vcl_string > to_return;

  bool bLoadComps;
  parameters()->get_value( "-bLoadComps", bLoadComps);

  if (bLoadComps){ //load 3 component images
    to_return.push_back( "image" );
    to_return.push_back( "image" );
    to_return.push_back( "image" );
  }
  else //single RGB image
    to_return.push_back( "image" );

  return to_return;
}


//: Provide a vector of output types
vcl_vector< vcl_string > dbdet_third_order_color_edge_detector_process::get_output_type()
{
  vcl_vector<vcl_string > to_return;
  to_return.push_back( "edge_map" );
  return to_return;
}


//: Execute the process
bool
dbdet_third_order_color_edge_detector_process::execute()
{
  //1) get the parameters
  bool bLoadComps, reduce_tokens;
  unsigned grad_op, conv_algo, col_conv_opt, parabola_fit, out_type;
  double sigma, thresh;
  int N=0;
  
  parameters()->get_value( "-bLoadComps", bLoadComps);
  parameters()->get_value( "-col_conv", col_conv_opt);
  parameters()->get_value( "-grad_op", grad_op);
  parameters()->get_value( "-conv_algo", conv_algo);
  parameters()->get_value( "-int_factor" , N );
  parameters()->get_value( "-sigma", sigma);
  parameters()->get_value( "-thresh", thresh);
  parameters()->get_value( "-parabola_fit", parabola_fit );
  parameters()->get_value( "-out_type", out_type );
  parameters()->get_value( "-breduce" , reduce_tokens );

  if (bLoadComps && input_data_.size() != 3 ){
    vcl_cout << "In dbdet_third_order_color_edge_detector_process::execute() - 3 images needed! \n";
    return false;
  }
  else if (input_data_.size() != 1 ){
    vcl_cout << "In dbdet_third_order_color_edge_detector_process::execute() - only one color image needed!\n";
    return false;
  }
  clear_output();

  vcl_cout << "Third Order color edge detection...";
  vcl_cout.flush();

  //2) get image(s) from the storage class
  vil_image_resource_sptr col_image_sptr, comp_image1_sptr, comp_image2_sptr, comp_image3_sptr;
  vil_image_view<vxl_byte> col_image;
  vil_image_view<float> comp1, comp2, comp3;

  if (bLoadComps){
    vidpro1_image_storage_sptr frame_image;
    frame_image.vertical_cast(input_data_[0][0]);
    comp_image1_sptr = frame_image->get_image();
    comp1 = comp_image1_sptr->get_view(0, comp_image1_sptr->ni(), 0, comp_image1_sptr->nj() );

    frame_image.vertical_cast(input_data_[0][1]);
    comp_image2_sptr = frame_image->get_image();
    comp2 = comp_image2_sptr->get_view(0, comp_image2_sptr->ni(), 0, comp_image2_sptr->nj() );

    frame_image.vertical_cast(input_data_[0][2]);
    comp_image3_sptr = frame_image->get_image();
    comp3 = comp_image3_sptr->get_view(0, comp_image3_sptr->ni(), 0, comp_image3_sptr->nj() );

    //make sure these images are one plane images
    if (comp1.nplanes() != 1 || comp1.nplanes() != 1 || comp1.nplanes() != 1){
      vcl_cout << "In dbdet_third_order_color_edge_detector_process::execute() - component images must be monochromatic! \n";
      return false;
    }
  }
  else {
    vidpro1_image_storage_sptr frame_image;
    frame_image.vertical_cast(input_data_[0][0]);
    col_image_sptr = frame_image->get_image();
    col_image = col_image_sptr->get_view(0, col_image_sptr->ni(), 0, col_image_sptr->nj() );

    //make sure these images are one plane images
    if (col_image.nplanes() != 3){
      vcl_cout << "In dbdet_third_order_color_edge_detector_process::execute() - image must be trichromatic! \n";
      return false;
    }
  }

  int padding=10;
  vil_image_view<vxl_byte> padded_img;
  padded_img.set_size(
      col_image.ni()+2*padding,col_image.nj()+2*padding,
      col_image.nplanes());
  vil_fill(padded_img, vxl_byte(0));

  vil_border_accessor<vil_image_view<vxl_byte> >
    accessor = vil_border_create_accessor(
        col_image,
        vil_border_create_geodesic(col_image));

  int j_max = (int)(padded_img.nj())-padding;
  int i_max = (int)(padded_img.ni())-padding;

  for (int p=0;p<padded_img.nplanes();++p)
  {
      for (int j = -padding ; j < j_max;++j)
      {
          for (int i=-padding;i < i_max;++i)
          {                          
              padded_img(i+padding,j+padding,p)=accessor(i,j,p);
          }
      }
  }
 

  //3) convert to the desired color space
  if (!bLoadComps){
    switch(col_conv_opt){
      case 0: //RGB color space
        vil_convert_cast(vil_plane(padded_img, 0), comp1);
        vil_convert_cast(vil_plane(padded_img, 1), comp2);
        vil_convert_cast(vil_plane(padded_img, 2), comp3);
        break;
      case 1: //IHS color space
        brip_vil_float_ops::convert_to_IHS(padded_img, comp1, comp2, comp3);
        break;
      case 2: //Lab color space
        convert_RGB_to_Lab(padded_img, comp1, comp2, comp3);
        break;
      case 3: //Luv color space
        convert_RGB_to_Luv(padded_img, comp1, comp2, comp3);
        break;
    }
  }

  //start the timer
  vul_timer t;

  dbdet_edgemap_sptr edge_map = dbdet_third_order_color(grad_op, conv_algo, N, sigma, thresh, parabola_fit, comp1, comp2, comp3, reduce_tokens);
  

  //create a new edgemap from the tokens collected from NMS
  dbdet_edgemap_sptr padded_edge_map = new dbdet_edgemap(col_image.ni(), 
                                                  col_image.nj());

  vcl_vector<dbdet_edgel*> padded_edges=edge_map->edgels;
  for ( unsigned int i=0; i < padded_edges.size() ; ++i)
  {

      vgl_point_2d<double> new_location;
      new_location.set(padded_edges[i]->pt.x()-(double)padding,
                       padded_edges[i]->pt.y()-(double)padding);
      padded_edges[i]->pt.set(new_location.x(),new_location.y());

      if ( (new_location.x() >= 0) && 
           (new_location.x() <= (double)col_image.ni()-1) && 
           (new_location.y() >= 0) &&
           (new_location.y() <= (double)col_image.nj()-1))
      {
          
          padded_edge_map->
              insert(new dbdet_edgel(padded_edges[i]->pt,
                                     padded_edges[i]->tangent,
                                     padded_edges[i]->strength,
                                     padded_edges[i]->deriv));
          
      } 
  }
//  edge_map->unref();
  edge_map=0;
  //double third_order_time = t.real();
  
  vcl_cout << "done!" << vcl_endl;
  
  //vcl_cout << "time taken for third_order: " << third_order_time << " msec" << vcl_endl;
  vcl_cout << "#edgels = " << padded_edge_map->num_edgels() << vcl_endl;

  // create the output storage class
  dbdet_edgemap_storage_sptr output_edgemap = dbdet_edgemap_storage_new();
  output_edgemap->set_edgemap(padded_edge_map);
  output_data_[0].push_back(output_edgemap);

  return true;
}

bool
dbdet_third_order_color_edge_detector_process::finish()
{
  return true;
}
