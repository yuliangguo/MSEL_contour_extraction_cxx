// This is brcv/seg/dbdet/algo/dbdet_subpix_convolution.h
#ifndef dbdet_subpix_convolution_h
#define dbdet_subpix_convolution_h
//:
//\file
//\brief A class to compute subpixel convolutions by using shifted operators
//\author Amir Tamrakar
//\date 09/09/06
//
//\verbatim
//  Modifications
//\endverbatim

#include <vil/algo/vil_correlate_1d.h>
#include <vil/algo/vil_correlate_2d.h>
#include <vil/vil_transpose.h>
#include <vgl/vgl_point_2d.h>

//: Convolve an image(src_im) with shifted kernels to compute the value at the subpixel positions
//  The subpixel values are computed at 1/2^N intervals in each dimension. e.g., N=0: regular
//  convolution, N=1: compute at mid-pixel locations, N=2: compute at quarter pixel locations
//
// I'm using vil_correlate_2d_at_pt to evaluate the dot_product at pixel/subpixel locations
// dest_im is resized to [2^N*src_im.ni() x 2^N*src_im.nj()]
//
template <class srcT, class destT, class kernelT, class accumT>
inline void dbdet_subpix_convolve_2d(const vil_image_view<srcT>& src_im,
                                     vil_image_view<destT>& dest_im,
                                     kernelT kernel, accumT ac, 
                                     int N)
{
  int two_pow_N = (int) vcl_pow(2.0, N);
  
  //the portion of the src_image that can be used for the convolution
  int src_ni = 1+src_im.ni()-kernel.ni(); assert(src_ni >= 0);
  int src_nj = 1+src_im.nj()-kernel.nj(); assert(src_nj >= 0);

  // filtered image will be 2^N time larger than the original
  int dest_ni = two_pow_N*src_im.ni(); assert(dest_ni >= 0);
  int dest_nj = two_pow_N*src_im.nj(); assert(dest_nj >= 0);

  //resize the dest image to this size
  dest_im.set_size(dest_ni,dest_nj,1);
  dest_im.fill(0);

  //determine the step sizes
  vcl_ptrdiff_t s_istep = src_im.istep(), s_jstep = src_im.jstep();
  vcl_ptrdiff_t s_pstep = src_im.planestep();
  vcl_ptrdiff_t d_istep = dest_im.istep(), d_jstep = dest_im.jstep();

  // filter with shifted versions of the original filter
  // we need 2^(2N) convolutions to completely compute it
  
  for (int shift_i=0; shift_i<two_pow_N; shift_i++){
    for (int shift_j=0; shift_j<two_pow_N; shift_j++){
      //instantiate a shifted kernel
      kernel.recompute_kernel(double(shift_i)/two_pow_N, double(shift_j)/two_pow_N);
      
      //invert kernel to use correlation code
      //kernelT inv_kernel = vil_flip_ud(vil_flip_lr(kernel));

      //***********************************************************
      //Now convolve the image with this kernel
      
      //Put the result in the right place on the destination image
      //depending on the current shifts

      //***********************************************************

      //Select the first pixel on the src image that is valid
      int khs = (kernel.ni()-1)/2;
      const srcT*  src_row  = src_im.top_left_ptr();

      //select the corresponding first pixel on the dest image
      destT* dest_row = dest_im.top_left_ptr()+(shift_i+two_pow_N*khs)+(shift_j+two_pow_N*khs)*d_jstep;

      for (int j=0;j<src_nj;++j,src_row+=s_jstep,dest_row+=two_pow_N*d_jstep)
      {
        const srcT* sp = src_row;
        destT* dp = dest_row;
        for (int i=0;i<src_ni;++i, sp += s_istep, dp += two_pow_N*d_istep){
          // Correlate at src(i,j)
          *dp = (destT)vil_correlate_2d_at_pt(sp, s_istep, s_jstep, s_pstep, kernel, ac);
        }
      }
    }
  }
}

//: Convolve the image with the given kernel at the given subpixel locations (pts) only
template <class srcT, class kernelT, class accumT>
inline void dbdet_subpix_convolve_2d(const vil_image_view<srcT>& src_im,
                                     const vcl_vector<vgl_point_2d<accumT> >& pts,
                                     vcl_vector<accumT>& res,
                                     kernelT kernel, accumT ac, 
                                     int /*N*/)
{
  res.resize(pts.size());
  int khs = (kernel.ni()-1)/2;

  //determine the step sizes for the src_image
  vcl_ptrdiff_t s_istep = src_im.istep(), s_jstep = src_im.jstep(), s_pstep = src_im.planestep();
  //for each subpixel point, we need to create a new kernel shifted to that location
  for (unsigned i=0; i<pts.size(); i++)
  {
    //determine closest integer coordinates of the current point
    int x = (int) vcl_floor(pts[i].x());
    int y = (int) vcl_floor(pts[i].y());

    //assert(x>=khs && y>=khs && x<src_im.ni()-khs && y<src_im.nj()-khs);

    if (x<khs || y<khs){
      res[i]=0;
      continue;
    }

    //determine the shift from it
    double dx = pts[i].x() - x;
    double dy = pts[i].y() - y;

    //instantiate a shifted kernel for this point
    kernel.recompute_kernel(dx, dy);
      
    //invert kernel to use correlation code
    //kernelT inv_kernel = vil_flip_ud(vil_flip_lr(kernel));

    //Select the first pixel on the src image that is valid
    
    const srcT*  sp  = src_im.top_left_ptr() + (x-khs)*s_istep + (y-khs)*s_jstep;
    
    // Correlate at src(x,y)
    res[i] = (accumT)vil_correlate_2d_at_pt(sp, s_istep, s_jstep, s_pstep, kernel, ac);
  }
}

//: Convolve the image with the given kernel at the given subpixel locations (pts) only at the given orientation
template <class srcT, class kernelT, class accumT>
inline void dbdet_subpix_convolve_2d(const vil_image_view<srcT>& src_im,
                                     const vcl_vector<vgl_point_2d<accumT> >& pts,
                                     const vcl_vector<accumT>& thetas,
                                     vcl_vector<accumT>& res,
                                     kernelT kernel, accumT ac, 
                                     int /*N*/)
{
  res.resize(pts.size());
  int khs = (kernel.ni()-1)/2;

  //determine the step sizes for the src_image
  vcl_ptrdiff_t s_istep = src_im.istep(), s_jstep = src_im.jstep(), s_pstep = src_im.planestep();
  //for each subpixel point, we need to create a new kernel shifted to that location
  for (unsigned i=0; i<pts.size(); i++)
  {
    //determine closest integer coordinates of the current point
    int x = (int) vcl_floor(pts[i].x());
    int y = (int) vcl_floor(pts[i].y());

    //assert(x>=khs && y>=khs && x<src_im.ni()-khs && y<src_im.nj()-khs);

    if (x<khs || y<khs){
      res[i]=0;
      continue;
    }

    //determine the shift from it
    double dx = pts[i].x() - x;
    double dy = pts[i].y() - y;

    //instantiate a shifted kernel for this point
    kernel.recompute_kernel(dx, dy, thetas[i]);
      
    //invert kernel to use correlation code
    //kernelT inv_kernel = vil_flip_ud(vil_flip_lr(kernel));

    //Select the first pixel on the src image that is valid
    
    const srcT*  sp  = src_im.top_left_ptr() + (x-khs)*s_istep + (y-khs)*s_jstep;
    
    // Correlate at src(x,y)
    res[i] = (accumT)vil_correlate_2d_at_pt(sp, s_istep, s_jstep, s_pstep, kernel, ac);
  }
}

//: Convolve an image(src_im) with shifted kernels to compute the value at the subpixel positions
//  The subpixel values are computed at 1/2^N intervals in each dimension. e.g., N=0: regular
//  convolution, N=1: compute at mid-pixel locations, N=2: compute at quarter pixel locations
//
//  To speed up the convolution, I'm performing 2 independent convolutions with 1-d kernels.
//  dest_im is resized to [2^N*src_im.ni() x 2^N*src_im.nj()]
//
template <class srcT, class destT, class kernelT, class accumT>
inline void dbdet_subpix_convolve_2d_sep(const vil_image_view<srcT>& src_im,
                                         vil_image_view<destT>& dest_im,
                                         kernelT kernel, accumT ac, 
                                         int N)
{
  int two_pow_N = (int) vcl_pow(2.0, N);
  
  //the portion of the src_image that can be used for the convolution
  int src_ni = src_im.ni(); assert(src_ni >= 0);
  int src_nj = src_im.nj(); assert(src_nj >= 0);

  // filtered image will be 2^N time larger than the original
  int dest_ni = two_pow_N*src_im.ni(); assert(dest_ni >= 0);
  int dest_nj = two_pow_N*src_im.nj(); assert(dest_nj >= 0);

  //resize the dest image to this size
  dest_im.set_size(dest_ni,dest_nj,1);
  dest_im.fill(0);

  //determine the step sizes
  vcl_ptrdiff_t s_istep = src_im.istep(), s_jstep = src_im.jstep();
  vcl_ptrdiff_t d_istep = dest_im.istep(), d_jstep = dest_im.jstep();

  // filter with shifted versions of the original filter
  // we need 2^(2N) convolutions to completely compute it
  
  for (int shift_i=0; shift_i<two_pow_N; shift_i++){
    for (int shift_j=0; shift_j<two_pow_N; shift_j++){
      
      // Generate filter: instantiate a shifted kernel
      kernel.recompute_kernel(double(shift_i)/two_pow_N, double(shift_j)/two_pow_N);
      
      //determine the half width of the kernel
      int khs = kernel.khs;

      // Apply 1D convolution along i direction using the K_x filter
      vil_image_view<destT> work_im;
      vil_correlate_1d(src_im, work_im, &kernel.K_x[khs], -khs, khs,
                      ac, vil_convolve_no_extend, vil_convolve_no_extend);

      // Apply 1D convolution along j direction by applying K_y filter to its transpose
      vil_image_view<destT> work_im2;
      work_im2.set_size(src_im.ni(), src_im.nj(), src_im.nplanes());

      vil_image_view<destT> work_im_t = vil_transpose(work_im);
      vil_image_view<destT> work_im2_t = vil_transpose(work_im2);

      vil_correlate_1d(work_im_t, work_im2_t,
                      &kernel.K_y[khs], -khs, khs,
                      ac, vil_convolve_no_extend, vil_convolve_no_extend);

      //compose this convolution, at the sub pixel location, into the larger dest image
      destT*  src_row  = work_im2.top_left_ptr();

      //select the corresponding first pixel on the dest image
      destT* dest_row = dest_im.top_left_ptr()+(shift_i)+(shift_j)*d_jstep;

      for (int j=0;j<src_nj;++j, src_row+=s_jstep, dest_row+=two_pow_N*d_jstep)
      {
        const destT* sp = src_row;
        destT* dp = dest_row;
        for (int i=0;i<src_ni;++i, sp += s_istep, dp += two_pow_N*d_istep){
          *dp = (destT) *sp;
        }
      }
    }
  }


}

////: correlate kernel[i] (i in [k_lo,k_hi]) with srcT in i-direction
//// On exit dest_im(i,j) = sum src(i+x,j)*kernel(x)  (x=k_lo..k_hi)
//// \note  This function does not reverse the kernel. If you want the
//// kernel reversed, use vil_convolve_1d instead.
//// \param kernel should point to tap 0.
//// \param dest_im will be resized to size of src_im.
//// \relates vil_image_view
//template <class srcT, class destT, class kernelT, class accumT>
//inline void dbdet_correlate_1d(const vil_image_view<srcT>& src_im,
//                             vil_image_view<destT>& dest_im,
//                             const kernelT* kernel,
//                             vcl_ptrdiff_t k_lo, vcl_ptrdiff_t k_hi,
//                             accumT ac,
//                             vil_convolve_boundary_option start_option,
//                             vil_convolve_boundary_option end_option)
//{
//  unsigned ni = src_im.ni();
//  unsigned nj = src_im.nj();
//  vcl_ptrdiff_t s_istep = src_im.istep(), s_jstep = src_im.jstep();
//
//  dest_im.set_size(ni,nj,src_im.nplanes());
//  vcl_ptrdiff_t d_istep = dest_im.istep(),d_jstep = dest_im.jstep();
//
//  for (unsigned int p=0;p<src_im.nplanes();++p)
//  {
//    // Select first row of p-th plane
//    const srcT*  src_row  = src_im.top_left_ptr()+p*src_im.planestep();
//    destT* dest_row = dest_im.top_left_ptr()+p*dest_im.planestep();
//
//    // Apply convolution to each row in turn
//    // First check if either istep is 1 for speed optimisation.
//    if (s_istep == 1)
//    {
//      if (d_istep == 1)
//        for (unsigned int j=0;j<nj;++j,src_row+=s_jstep,dest_row+=d_jstep)
//          vil_correlate_1d(src_row,ni,1,  dest_row,1,
//                           kernel,k_lo,k_hi,ac,start_option,end_option);
//      else
//        for (unsigned int j=0;j<nj;++j,src_row+=s_jstep,dest_row+=d_jstep)
//          vil_correlate_1d(src_row,ni,1,  dest_row,d_istep,
//                           kernel,k_lo,k_hi,ac,start_option,end_option);
//    }
//    else
//    {
//      if (d_istep == 1)
//        for (unsigned int j=0;j<nj;++j,src_row+=s_jstep,dest_row+=d_jstep)
//          vil_correlate_1d(src_row,ni,s_istep,  dest_row,1,
//                           kernel,k_lo,k_hi,ac,start_option,end_option);
//      else
//        for (unsigned int j=0;j<nj;++j,src_row+=s_jstep,dest_row+=d_jstep)
//          vil_correlate_1d(src_row,ni,s_istep,  dest_row,d_istep,
//                           kernel,k_lo,k_hi,ac,start_option,end_option);
//    }
//  }
//}

#endif // dbdet_subpix_convolution_h
