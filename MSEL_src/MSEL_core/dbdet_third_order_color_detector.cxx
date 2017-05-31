#include "dbdet_third_order_color_detector.h"

#include <vul/vul_timer.h>
#include "dbdet_edgemap.h"

#include "dbdet_gaussian_kernel.h"
#include "dbdet_interp_kernel.h"
#include "dbdet_nms.h"
#include "dbdet_subpix_convolution.h"
#include "dbdet_edgemap_sptr.h"
#include "dbdet_sel_utils.h"
#include "dbdet_edge_appearance_util.h"

#include <vil/algo/vil_sobel_1x3.h>
#include <vil/algo/vil_sobel_3x3.h>
#include <vil/algo/vil_convolve_2d.h>
#include <vil/vil_resample_bilin.h>
#include <brip/brip_vil_float_ops.h>


dbdet_edgemap_sptr dbdet_third_order_color(unsigned grad_op, unsigned conv_algo,
                                           int N, double sigma, double thresh, unsigned parabola_fit,
                                           vil_image_view<float>& comp1, 
                                           vil_image_view<float>& comp2, 
                                           vil_image_view<float>& comp3,
                                           bool reduce_tokens)
{
  //start the timer
  vul_timer t;

  //4) compute image gradients of each of the components
  vil_image_view<float> f1_dx, f1_dy, f2_dx, f2_dy, f3_dx, f3_dy;
  int scale=1;

  switch (grad_op)
  {
    case 0: //Gaussian
    {  
      scale = (int) vcl_pow(2.0, N);

      //compute gradients
      if (conv_algo==0){
        dbdet_subpix_convolve_2d(comp1, f1_dx, dbdet_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d(comp1, f1_dy, dbdet_Gy_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d(comp2, f2_dx, dbdet_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d(comp2, f2_dy, dbdet_Gy_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d(comp3, f3_dx, dbdet_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d(comp3, f3_dy, dbdet_Gy_kernel(sigma), float(), N);
      }
      else {
        dbdet_subpix_convolve_2d_sep(comp1, f1_dx, dbdet_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d_sep(comp1, f1_dy, dbdet_Gy_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d_sep(comp2, f2_dx, dbdet_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d_sep(comp2, f2_dy, dbdet_Gy_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d_sep(comp3, f3_dx, dbdet_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d_sep(comp3, f3_dy, dbdet_Gy_kernel(sigma), float(), N);
      }
      break;
    }
    case 1: //h0-operator
    {
      scale = (int) vcl_pow(2.0, N);

      //compute gradients  
      if (conv_algo==0){
        dbdet_subpix_convolve_2d(comp1, f1_dx, dbdet_h0_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d(comp1, f1_dy, dbdet_h0_Gy_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d(comp2, f2_dx, dbdet_h0_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d(comp2, f2_dy, dbdet_h0_Gy_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d(comp3, f3_dx, dbdet_h0_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d(comp3, f3_dy, dbdet_h0_Gy_kernel(sigma), float(), N);
      }
      else {
        dbdet_subpix_convolve_2d_sep(comp1, f1_dx, dbdet_h0_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d_sep(comp1, f1_dy, dbdet_h0_Gy_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d_sep(comp2, f2_dx, dbdet_h0_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d_sep(comp2, f2_dy, dbdet_h0_Gy_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d_sep(comp3, f3_dx, dbdet_h0_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d_sep(comp3, f3_dy, dbdet_h0_Gy_kernel(sigma), float(), N);
      }

      break;
    }
    case 2:  //h1-operator
    {
      scale = (int) vcl_pow(2.0, N);

      //compute gradients
      if (conv_algo==0){
        dbdet_subpix_convolve_2d(comp1, f1_dx, dbdet_h1_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d(comp1, f1_dy, dbdet_h1_Gy_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d(comp2, f2_dx, dbdet_h1_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d(comp2, f2_dy, dbdet_h1_Gy_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d(comp3, f3_dx, dbdet_h1_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d(comp3, f3_dy, dbdet_h1_Gy_kernel(sigma), float(), N);
      }
      else {
        dbdet_subpix_convolve_2d_sep(comp1, f1_dx, dbdet_h1_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d_sep(comp1, f1_dy, dbdet_h1_Gy_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d_sep(comp2, f2_dx, dbdet_h1_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d_sep(comp2, f2_dy, dbdet_h1_Gy_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d_sep(comp3, f3_dx, dbdet_h1_Gx_kernel(sigma), float(), N);
        dbdet_subpix_convolve_2d_sep(comp3, f3_dy, dbdet_h1_Gy_kernel(sigma), float(), N);
      }
      break;
    }
  }
  
  //5) compute the squared norm of the vector-gradient
  vil_image_view<double> grad_mag, nu1, nu2; //eigenvalue and eigenvector
  grad_mag.set_size(f1_dx.ni(), f1_dx.nj());
  nu1.set_size(f1_dx.ni(), f1_dx.nj());
  nu2.set_size(f1_dx.ni(), f1_dx.nj());

  //get the pointers to the memory chunks
  float *f1x  =  f1_dx.top_left_ptr();
  float *f1y  =  f1_dy.top_left_ptr();
  float *f2x  =  f2_dx.top_left_ptr();
  float *f2y  =  f2_dy.top_left_ptr();
  float *f3x  =  f3_dx.top_left_ptr();
  float *f3y  =  f3_dy.top_left_ptr();
  double *g_mag  =  grad_mag.top_left_ptr();
  double *n1  =  nu1.top_left_ptr();
  double *n2  =  nu2.top_left_ptr();

  //compute the squared norm of gradient
  for(unsigned long i=0; i<grad_mag.size(); i++){
    double A = f1x[i]*f1x[i]+f2x[i]*f2x[i]+f3x[i]*f3x[i];
    double B = f1x[i]*f1y[i]+f2x[i]*f2y[i]+f3x[i]*f3y[i];
    double C = f1y[i]*f1y[i]+f2y[i]*f2y[i]+f3y[i]*f3y[i];

    double l = (A+C+ vcl_sqrt((A-C)*(A-C) + 4*B*B))/2;

    if (vcl_fabs(B)>1e-2){
      n1[i] = B/vcl_sqrt(B*B+(l-A)*(l-A));
      n2[i] = (l-A)/vcl_sqrt(B*B+(l-A)*(l-A));
    }
    else {
      n1[i] = (l-C)/vcl_sqrt(B*B+(l-C)*(l-C));
      n2[i] = B/vcl_sqrt(B*B+(l-C)*(l-C));
    }

    g_mag[i] = vcl_sqrt(l/3); //take the square root of the squared norm
  }


  double conv_time = t.real();  
  t.mark(); //reset timer

  //Now call the nms code to get the subpixel edge tokens
  vcl_vector<vgl_point_2d<double> > edge_locations;
  vcl_vector<double> orientation, mag, d2f;

  dbdet_nms NMS(dbdet_nms_params(thresh, (dbdet_nms_params::PFIT_TYPE)parabola_fit), nu1, nu2, grad_mag);
  NMS.apply(true, edge_locations, orientation, mag, d2f);

  double nms_time = t.real();
  t.mark(); //reset timer

  //convert to the original image scale coordinates
  for (unsigned i=0; i<edge_locations.size(); i++)
    edge_locations[i].set(edge_locations[i].x()/scale, edge_locations[i].y()/scale);

  //for each edge, compute all the gradients to compute the new orientation
  vcl_vector<double> If1x, If1y, If1xx, If1xy, If1yy, If1xxy, If1xyy, If1xxx, If1yyy;
  vcl_vector<double> If2x, If2y, If2xx, If2xy, If2yy, If2xxy, If2xyy, If2xxx, If2yyy;
  vcl_vector<double> If3x, If3y, If3xx, If3xy, If3yy, If3xxy, If3xyy, If3xxx, If3yyy;

  switch (grad_op)
  {
    case 0: //Interpolated Gaussian
    { 
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1x,   dbdet_Gx_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1y,   dbdet_Gy_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1xx,  dbdet_Gxx_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1xy,  dbdet_Gxy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1yy,  dbdet_Gyy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1xxx, dbdet_Gxxx_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1xxy, dbdet_Gxxy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1xyy, dbdet_Gxyy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1yyy, dbdet_Gyyy_kernel(sigma), double(), N);

      dbdet_subpix_convolve_2d(comp2, edge_locations, If2x,   dbdet_Gx_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2y,   dbdet_Gy_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2xx,  dbdet_Gxx_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2xy,  dbdet_Gxy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2yy,  dbdet_Gyy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2xxx, dbdet_Gxxx_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2xxy, dbdet_Gxxy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2xyy, dbdet_Gxyy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2yyy, dbdet_Gyyy_kernel(sigma), double(), N);

      dbdet_subpix_convolve_2d(comp3, edge_locations, If3x,   dbdet_Gx_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3y,   dbdet_Gy_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3xx,  dbdet_Gxx_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3xy,  dbdet_Gxy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3yy,  dbdet_Gyy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3xxx, dbdet_Gxxx_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3xxy, dbdet_Gxxy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3xyy, dbdet_Gxyy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3yyy, dbdet_Gyyy_kernel(sigma), double(), N);

      break;
    }
    case 1: //h0-operator
    {
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1x,   dbdet_h0_Gx_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1y,   dbdet_h0_Gy_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1xx,  dbdet_h0_Gxx_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1xy,  dbdet_h0_Gxy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1yy,  dbdet_h0_Gyy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1xxx, dbdet_h0_Gxxx_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1xxy, dbdet_h0_Gxxy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1xyy, dbdet_h0_Gxyy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1yyy, dbdet_h0_Gyyy_kernel(sigma), double(), N);

      dbdet_subpix_convolve_2d(comp2, edge_locations, If2x,   dbdet_h0_Gx_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2y,   dbdet_h0_Gy_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2xx,  dbdet_h0_Gxx_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2xy,  dbdet_h0_Gxy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2yy,  dbdet_h0_Gyy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2xxx, dbdet_h0_Gxxx_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2xxy, dbdet_h0_Gxxy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2xyy, dbdet_h0_Gxyy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2yyy, dbdet_h0_Gyyy_kernel(sigma), double(), N);

      dbdet_subpix_convolve_2d(comp3, edge_locations, If3x,   dbdet_h0_Gx_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3y,   dbdet_h0_Gy_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3xx,  dbdet_h0_Gxx_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3xy,  dbdet_h0_Gxy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3yy,  dbdet_h0_Gyy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3xxx, dbdet_h0_Gxxx_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3xxy, dbdet_h0_Gxxy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3xyy, dbdet_h0_Gxyy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3yyy, dbdet_h0_Gyyy_kernel(sigma), double(), N);

      break;
    }
    case 2:  //h1-operator
    {
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1x,   dbdet_h1_Gx_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1y,   dbdet_h1_Gy_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1xx,  dbdet_h1_Gxx_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1xy,  dbdet_h1_Gxy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1yy,  dbdet_h1_Gyy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1xxx, dbdet_h1_Gxxx_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1xxy, dbdet_h1_Gxxy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1xyy, dbdet_h1_Gxyy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp1, edge_locations, If1yyy, dbdet_h1_Gyyy_kernel(sigma), double(), N);

      dbdet_subpix_convolve_2d(comp2, edge_locations, If2x,   dbdet_h1_Gx_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2y,   dbdet_h1_Gy_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2xx,  dbdet_h1_Gxx_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2xy,  dbdet_h1_Gxy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2yy,  dbdet_h1_Gyy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2xxx, dbdet_h1_Gxxx_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2xxy, dbdet_h1_Gxxy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2xyy, dbdet_h1_Gxyy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp2, edge_locations, If2yyy, dbdet_h1_Gyyy_kernel(sigma), double(), N);

      dbdet_subpix_convolve_2d(comp3, edge_locations, If3x,   dbdet_h1_Gx_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3y,   dbdet_h1_Gy_kernel(sigma),   double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3xx,  dbdet_h1_Gxx_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3xy,  dbdet_h1_Gxy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3yy,  dbdet_h1_Gyy_kernel(sigma),  double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3xxx, dbdet_h1_Gxxx_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3xxy, dbdet_h1_Gxxy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3xyy, dbdet_h1_Gxyy_kernel(sigma), double(), N);
      dbdet_subpix_convolve_2d(comp3, edge_locations, If3yyy, dbdet_h1_Gyyy_kernel(sigma), double(), N);

      break;
    }
  }
      
  //Now, compute and update each edge with its new orientation
  vcl_vector<double> edge_orientations(edge_locations.size());
  
  for (unsigned i=0; i<edge_locations.size();i++)
  {
    double A   = If1x[i]*If1x[i] + If2x[i]*If2x[i] + If3x[i]*If3x[i];
    double Ax  = 2*If1x[i]*If1xx[i] + 2*If2x[i]*If2xx[i] + 2*If3x[i]*If3xx[i];
    double Ay  = 2*If1x[i]*If1xy[i] + 2*If2x[i]*If2xy[i] + 2*If3x[i]*If3xy[i];
    double Axx = 2*If1x[i]*If1xxx[i] + 2*If1xx[i]*If1xx[i] + 2*If2x[i]*If2xxx[i] + 2*If2xx[i]*If2xx[i] + 2*If3x[i]*If3xxx[i] + 2*If3xx[i]*If3xx[i];
    double Axy = 2*If1x[i]*If1xxy[i] + 2*If1xy[i]*If1xx[i] + 2*If2x[i]*If2xxy[i] + 2*If2xy[i]*If2xx[i] + 2*If3x[i]*If3xxy[i] + 2*If3xy[i]*If3xx[i];
    double Ayy = 2*If1x[i]*If1xyy[i] + 2*If1xy[i]*If1xy[i] + 2*If2x[i]*If2xyy[i] + 2*If2xy[i]*If2xy[i] + 2*If3x[i]*If3xyy[i] + 2*If3xy[i]*If3xy[i];

    double B   = If1x[i]*If1y[i] + If2x[i]*If2y[i] + If3x[i]*If3y[i];
    double Bx  = If1x[i]*If1xy[i] + If1y[i]*If1xx[i] + If2x[i]*If2xy[i] + If2y[i]*If2xx[i] + If3x[i]*If3xy[i] + If3y[i]*If3xx[i];
    double By  = If1x[i]*If1yy[i] + If1y[i]*If1xy[i] + If2x[i]*If2yy[i] + If2y[i]*If2xy[i] + If3x[i]*If3yy[i] + If3y[i]*If3xy[i];
    double Bxx = If1x[i]*If1xxy[i] + If1xx[i]*If1xy[i] + If1y[i]*If1xxx[i] + If1xy[i]*If1xx[i] + If2x[i]*If2xxy[i] + If2xx[i]*If2xy[i] + If2y[i]*If2xxx[i] + If2xy[i]*If2xx[i] + If3x[i]*If3xxy[i] + If3xx[i]*If3xy[i] + If3y[i]*If3xxx[i] + If3xy[i]*If3xx[i];
    double Bxy = If1x[i]*If1xyy[i] + If1xy[i]*If1xy[i] + If1y[i]*If1xxy[i] + If1yy[i]*If1xx[i] + If2x[i]*If2xyy[i] + If2xy[i]*If2xy[i] + If2y[i]*If2xxy[i] + If2yy[i]*If2xx[i] + If3x[i]*If3xyy[i] + If3xy[i]*If3xy[i] + If3y[i]*If3xxy[i] + If3yy[i]*If3xx[i];
    double Byy = If1x[i]*If1yyy[i] + If1xy[i]*If1yy[i] + If1y[i]*If1xyy[i] + If1yy[i]*If1xy[i] + If2x[i]*If2yyy[i] + If2xy[i]*If2yy[i] + If2y[i]*If2xyy[i] + If2yy[i]*If2xy[i] + If3x[i]*If3yyy[i] + If3xy[i]*If3yy[i] + If3y[i]*If3xyy[i] + If3yy[i]*If3xy[i];

    double C   = If1y[i]*If1y[i] + If2y[i]*If2y[i] + If3y[i]*If3y[i];
    double Cx  = 2*If1y[i]*If1xy[i] + 2*If2y[i]*If2xy[i] + 2*If3y[i]*If3xy[i];
    double Cy  = 2*If1y[i]*If1yy[i] + 2*If2y[i]*If2yy[i] + 2*If3y[i]*If3yy[i];
    double Cxx = 2*If1y[i]*If1xxy[i] + 2*If1xy[i]*If1xy[i] + 2*If2y[i]*If2xxy[i] + 2*If2xy[i]*If2xy[i] + 2*If3y[i]*If3xxy[i] + 2*If3xy[i]*If3xy[i];
    double Cxy = 2*If1y[i]*If1xyy[i] + 2*If1yy[i]*If1xy[i] + 2*If2y[i]*If2xyy[i] + 2*If2yy[i]*If2xy[i] + 2*If3y[i]*If3xyy[i] + 2*If3yy[i]*If3xy[i];
    double Cyy = 2*If1y[i]*If1yyy[i] + 2*If1yy[i]*If1yy[i] + 2*If2y[i]*If2yyy[i] + 2*If2yy[i]*If2yy[i] + 2*If3y[i]*If3yyy[i] + 2*If3yy[i]*If3yy[i];

    double l = ((C+A) + vcl_sqrt((A-C)*(A-C) + 4*B*B))/2.0;
    double e = (2*l-A-C);

    double lx = ( (l-C)*Ax + (l-A)*Cx + 2*B*Bx)/e;
    double ly = ( (l-C)*Ay + (l-A)*Cy + 2*B*By)/e;
    double lxx = ( (lx-Cx)*Ax + (l-C)*Axx + (lx-Ax)*Cx + (l-A)*Cxx + 2*Bx*Bx + 2*B*Bxx - lx*(2*lx-Cx-Ax))/e;
    double lxy = ( (ly-Cy)*Ax + (l-C)*Axy + (ly-Ay)*Cx + (l-A)*Cxy + 2*Bx*By + 2*B*Bxy - lx*(2*ly-Cy-Ay))/e;
    double lyy = ( (ly-Cy)*Ay + (l-C)*Ayy + (ly-Ay)*Cy + (l-A)*Cyy + 2*By*By + 2*B*Byy - ly*(2*ly-Cy-Ay))/e;

    /********************************************************************************
    // This computation is noisy
    double n1 = vcl_sqrt((1.0 + (A-C)/d)/2.0);
    double n2 = vnl_math_sgn(B)*vcl_sqrt((1.0 - (A-C)/d)/2.0);

    // when B is zero, these derivatives need to be corrected to the limiting value
    double n1x, n1y, n2x, n2y;
    if (vcl_fabs(B)>1e-2){
      n1x = ( 2*(Ax-Cx)/d - (A-C)*(2*(C-A)*(Cx-Ax) + 8*B*Bx)/(d*d*d))/2/n1;
      n1y = ( 2*(Ay-Cy)/d - (A-C)*(2*(C-A)*(Cy-Ay) + 8*B*By)/(d*d*d))/2/n1;
      n2x = (-2*(Ax-Cx)/d + (A-C)*(2*(C-A)*(Cx-Ax) + 8*B*Bx)/(d*d*d))/2/n2;
      n2y = (-2*(Ay-Cy)/d + (A-C)*(2*(C-A)*(Cy-Ay) + 8*B*By)/(d*d*d))/2/n2;
    }
    else {
      n1x = 0;
      n1y = 0;
      n2x = 0;
      n2y = 0; 
    }
    ********************************************************************************/

    double n1, n2, n1x, n1y, n2x, n2y;
    // when B is zero, these derivatives need to be fixed
    if (vcl_fabs(B)>1e-2){
      double f = vcl_sqrt(B*B+(l-A)*(l-A));

      n1 = B/f;
      n2 = (l-A)/f;

      n1x = Bx/f - B*(B*Bx + (l-A)*(lx-Ax))/(f*f*f);
      n1y = By/f - B*(B*By + (l-A)*(ly-Ay))/(f*f*f);
      n2x = (lx-Ax)/f - (l-A)*(B*Bx + (l-A)*(lx-Ax))/(f*f*f);
      n2y = (ly-Ay)/f - (l-A)*(B*By + (l-A)*(ly-Ay))/(f*f*f);
    }
    else {
      double f = vcl_sqrt(B*B+(l-C)*(l-C));

      n1 = (l-C)/f;
      n2 = B/f;

      n1x = (lx-Cx)/f - (l-C)*(B*Bx + (l-C)*(lx-Cx))/(f*f*f);
      n1y = (ly-Cy)/f - (l-C)*(B*By + (l-C)*(ly-Cy))/(f*f*f); 
      n2x = Bx/f - B*(B*Bx + (l-C)*(lx-Cx))/(f*f*f);
      n2y = By/f - B*(B*By + (l-C)*(ly-Cy))/(f*f*f);
    }

    //compute Fx and Fy
    double Fx = lx*n1x + n1*lxx + n2x*ly + n2*lxy;
    double Fy = lx*n1y + n1*lxy + n2y*ly + n2*lyy;

    //save new orientation
    edge_orientations[i] = dbdet_angle0To2Pi(vcl_atan2(Fx, -Fy));
  }

  double third_order_time = t.real();

  //report timings
  vcl_cout << vcl_endl;
  vcl_cout << "time taken for conv: " << conv_time << " msec" << vcl_endl;
  vcl_cout << "time taken for nms: " << nms_time << " msec" << vcl_endl;
  vcl_cout << "time taken for color third-order: " << third_order_time << " msec" << vcl_endl;
 

  //create a new edgemap from the tokens collected from NMS
  dbdet_edgemap_sptr edge_map = new dbdet_edgemap(comp1.ni(), comp1.nj());

  for (unsigned i=0; i<edge_locations.size(); i++){
    if (reduce_tokens){
      //only insert one edgel per grid position
      int xx = dbdet_round(edge_locations[i].x());
      int yy = dbdet_round(edge_locations[i].y());

      if (edge_map->edge_cells(yy, xx).size()>0)
        continue;

      edge_map->insert(new dbdet_edgel(edge_locations[i], edge_orientations[i], mag[i], d2f[i]));
    }
    else { //insert all of them
      edge_map->insert(new dbdet_edgel(edge_locations[i], edge_orientations[i], mag[i], d2f[i]));
    }
  }

  //add intensity appearance to the edges
  dbdet_add_color_app(comp1, comp2, comp3, edge_map, sigma, 1); //opt: 0: original , 1: smoothed, 2: Half gaussian

  return edge_map;
}

