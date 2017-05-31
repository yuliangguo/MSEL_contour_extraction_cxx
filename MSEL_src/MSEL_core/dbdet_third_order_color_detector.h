// This is brcv/seg/dbdet/algo/dbdet_third_order_color_detector.h
#ifndef dbdet_third_order_color_detector_h
#define dbdet_third_order_color_detector_h
//:
//\file
//\brief Utility functions for dbdet_third_order_color_detector
//\author Ozge Can Ozcanli adapted from Amir Tamrakar's code 
//\date 03/08/07 
//

#include "dbdet_edgemap_sptr.h"

#include <vil/vil_image_view.h>


dbdet_edgemap_sptr dbdet_third_order_color(unsigned grad_op, unsigned conv_algo,
                                           int N, double sigma, double thres, unsigned parabola_fit,
                                           vil_image_view<float>& comp1, 
                                           vil_image_view<float>& comp2, 
                                           vil_image_view<float>& comp3,
                                           bool reduce_tokens=false);

#endif // dbdet_third_order_color_detector_h

