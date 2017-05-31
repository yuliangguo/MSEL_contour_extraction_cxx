// This is brcv/seg/dbdet/algo/dbdet_edge_appearance_util.cxx

#include "dbdet_edge_appearance_util.h"

#include <vil/vil_image_resource.h>
#include <vil/vil_new.h>
#include <vil/vil_image_view.h>
#include <vil/vil_convert.h>
#include <vil/vil_bilin_interp.h>

#include "dbdet_gaussian_kernel.h"
#include "dbdet_interp_kernel.h"
#include "dbdet_nms.h"
#include "dbdet_subpix_convolution.h"
#include "dbdet_edgel.h"
#include "dbdet_edgemap.h"
#include "dbdet_edgemap_sptr.h"
#include "dbdet_sel_utils.h"

//: add intensity appearance to the edges
void dbdet_add_intensity_app(vil_image_view<vxl_byte>& greyscale_view, dbdet_edgemap_sptr edgemap, double sigma, int opt)
{
  switch(opt)
  {
    case 0 : //get intensity from the original image
    {
      for (unsigned i=0; i<edgemap->edgels.size(); i++)
      {
        dbdet_edgel* e = edgemap->edgels[i];

        //get the point to read off the intensity from
        vgl_point_2d<double> ptL(e->pt.x() - 1.5*sigma*vcl_cos(e->tangent-vnl_math::pi_over_2), 
                                 e->pt.y() - 1.5*sigma*vcl_sin(e->tangent-vnl_math::pi_over_2));
        vgl_point_2d<double> ptR(e->pt.x() + 1.5*sigma*vcl_cos(e->tangent-vnl_math::pi_over_2), 
                                 e->pt.y() + 1.5*sigma*vcl_sin(e->tangent-vnl_math::pi_over_2));

        // get the appearance data
        dbdet_appearance* lapp = new dbdet_intensity(vil_bilin_interp(greyscale_view, ptL.x(), ptL.y()));
        dbdet_appearance* rapp = new dbdet_intensity(vil_bilin_interp(greyscale_view, ptR.x(), ptR.y()));

        //add it to the edgel
        e->left_app = lapp;
        e->right_app = rapp;
      }
      break;
    }
    case 1: //get intensity from the smoothed image
    {
      //compute a Gaussian smoothed image from which to compute the intensities
      vil_image_view<double> sm_img;
      dbdet_subpix_convolve_2d(greyscale_view, sm_img, dbdet_G_kernel(sigma), double(), 0);

      for (unsigned i=0; i<edgemap->edgels.size(); i++)
      {
        dbdet_edgel* e = edgemap->edgels[i];

        //get the point to read off the intensity from
        vgl_point_2d<double> ptL(e->pt.x() - 1.5*sigma*vcl_cos(e->tangent-vnl_math::pi_over_2), 
                                 e->pt.y() - 1.5*sigma*vcl_sin(e->tangent-vnl_math::pi_over_2));
        vgl_point_2d<double> ptR(e->pt.x() + 1.5*sigma*vcl_cos(e->tangent-vnl_math::pi_over_2), 
                                 e->pt.y() + 1.5*sigma*vcl_sin(e->tangent-vnl_math::pi_over_2));

        // get the appearance data
        dbdet_appearance* lapp = new dbdet_intensity(vil_bilin_interp(sm_img, ptL.x(), ptL.y()));
        dbdet_appearance* rapp = new dbdet_intensity(vil_bilin_interp(sm_img, ptR.x(), ptR.y()));

        //add it to the edgel
        e->left_app = lapp;
        e->right_app = rapp;
      }
      break;
    }
    case 2: //get intensity from half Gaussian kernels
    {
      //compute the intensities on either side of the edge smoothed by a half gaussian
      vcl_vector<double> IL, IR;
      vcl_vector<vgl_point_2d<double> > pts(edgemap->edgels.size());
      vcl_vector<double> thetas(edgemap->edgels.size());
      
      for (unsigned i=0; i<edgemap->edgels.size(); i++)
      {
        dbdet_edgel* e = edgemap->edgels[i];
        pts[i] = e->pt;
        thetas[i] = e->tangent;
      }

      //compute the appearance values at these points
      dbdet_subpix_convolve_2d(greyscale_view, pts, thetas, IL, dbdet_G_Lhalf_kernel(2*sigma), double(), 1);
      dbdet_subpix_convolve_2d(greyscale_view, pts, thetas, IR, dbdet_G_Rhalf_kernel(2*sigma), double(), 1);

      for (unsigned i=0; i<edgemap->edgels.size(); i++)
      {
        dbdet_edgel* e = edgemap->edgels[i];

        //get the appearance data
        dbdet_appearance* lapp = new dbdet_intensity(IL[i]);
        dbdet_appearance* rapp = new dbdet_intensity(IR[i]);

        //add it to the edgel
        e->left_app = lapp;
        e->right_app = rapp;
      }
    }
    break;
  }
}

//: add color appearance to the edges
void dbdet_add_color_app(vil_image_view<float>& comp1, 
                         vil_image_view<float>& comp2, 
                         vil_image_view<float>& comp3,
                         dbdet_edgemap_sptr edgemap, double sigma, int opt)
{
  switch(opt)
  {
    case 0 : //get color from the original image
    {
      for (unsigned i=0; i<edgemap->edgels.size(); i++)
      {
        dbdet_edgel* e = edgemap->edgels[i];

        //get the point to read off the intensity from
        vgl_point_2d<double> ptL(e->pt.x() - 1.5*sigma*vcl_cos(e->tangent-vnl_math::pi_over_2), 
                                 e->pt.y() - 1.5*sigma*vcl_sin(e->tangent-vnl_math::pi_over_2));
        vgl_point_2d<double> ptR(e->pt.x() + 1.5*sigma*vcl_cos(e->tangent-vnl_math::pi_over_2), 
                                 e->pt.y() + 1.5*sigma*vcl_sin(e->tangent-vnl_math::pi_over_2));

        // get the appearance data
        dbdet_appearance* lapp = new dbdet_color(vil_bilin_interp(comp1, ptL.x(), ptL.y()),
                                                 vil_bilin_interp(comp2, ptL.x(), ptL.y()),
                                                 vil_bilin_interp(comp3, ptL.x(), ptL.y()));
        dbdet_appearance* rapp = new dbdet_color(vil_bilin_interp(comp1, ptR.x(), ptR.y()),
                                                 vil_bilin_interp(comp2, ptR.x(), ptR.y()),
                                                 vil_bilin_interp(comp3, ptR.x(), ptR.y()));

        //add it to the edgel
        e->left_app = lapp;
        e->right_app = rapp;
      }
      break;
    }
    case 1: //get color from the smoothed image
    {
      //compute a Gaussian smoothed image from which to compute the intensities
      vil_image_view<double> sm_comp1, sm_comp2, sm_comp3;
      dbdet_subpix_convolve_2d(comp1, sm_comp1, dbdet_G_kernel(sigma), double(), 0);
      dbdet_subpix_convolve_2d(comp2, sm_comp2, dbdet_G_kernel(sigma), double(), 0);
      dbdet_subpix_convolve_2d(comp3, sm_comp3, dbdet_G_kernel(sigma), double(), 0);

      for (unsigned i=0; i<edgemap->edgels.size(); i++)
      {
        dbdet_edgel* e = edgemap->edgels[i];

        //get the point to read off the intensity from
        vgl_point_2d<double> ptL(e->pt.x() - 1.5*sigma*vcl_cos(e->tangent-vnl_math::pi_over_2), 
                                 e->pt.y() - 1.5*sigma*vcl_sin(e->tangent-vnl_math::pi_over_2));
        vgl_point_2d<double> ptR(e->pt.x() + 1.5*sigma*vcl_cos(e->tangent-vnl_math::pi_over_2), 
                                 e->pt.y() + 1.5*sigma*vcl_sin(e->tangent-vnl_math::pi_over_2));

        // get the appearance data
        dbdet_appearance* lapp = new dbdet_color(vil_bilin_interp(sm_comp1, ptL.x(), ptL.y()),
                                                 vil_bilin_interp(sm_comp2, ptL.x(), ptL.y()),
                                                 vil_bilin_interp(sm_comp3, ptL.x(), ptL.y()));
        dbdet_appearance* rapp = new dbdet_color(vil_bilin_interp(sm_comp1, ptR.x(), ptR.y()),
                                                 vil_bilin_interp(sm_comp2, ptR.x(), ptR.y()),
                                                 vil_bilin_interp(sm_comp3, ptR.x(), ptR.y()));

        //add it to the edgel
        e->left_app = lapp;
        e->right_app = rapp;
      }
      break;
    }
  }
}

