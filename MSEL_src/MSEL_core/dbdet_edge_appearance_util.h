// This is brcv/seg/dbdet/algo/dbdet_edge_appearance_util.h
#ifndef dbdet_edge_appearance_util_h
#define dbdet_edge_appearance_util_h
//:
//\file
//\brief Functions to add appearance properties to the detected edges
//\author Amir Tamrakar
//\date 11/25/07
//
//\verbatim
//  Modifications
//\endverbatim

#include <vil/vil_image_view.h>
#include "dbdet_edgemap_sptr.h"

//: add intensity appearance to the edges
void dbdet_add_intensity_app(vil_image_view<vxl_byte>& greyscale_view, dbdet_edgemap_sptr edgemap, double sigma=1.0, int opt=0);

//: add color appearance to the edges
void dbdet_add_color_app(vil_image_view<float>& comp1, 
                         vil_image_view<float>& comp2, 
                         vil_image_view<float>& comp3,
                         dbdet_edgemap_sptr edgemap, double sigma=1.0, int opt=0);

#endif //dbdet_edge_appearance_util_h
