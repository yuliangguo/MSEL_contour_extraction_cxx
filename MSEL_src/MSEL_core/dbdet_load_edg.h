// This is brcv/seg/dbdet/algo/dbdet_load_edg.h
#ifndef dbdet_load_edg_h
#define dbdet_load_edg_h
//:
//\file
//\brief Utility functions for edge input
//\author Ricardo Fabbri (rfabbri), Brown University  (rfabbri.github.io)
//\date 01/11/07 15:56:05 EST
//
//\verbatim
//Modifications
//
//  Sep 29 2009  Ricardo Fabbri   Added support for writing gzip-compressed edg
//  Sep 16 2009  Ricardo Fabbri   Added support for reading gzip-compressed edg
//
//               Amir Tamrakar    Separated the load edge functions into two
//                                (for vsol and edgemap)
//
//               Amir Tamrakar    Modified to save and load grid locations as
//                                the integer coordinates.
//
// \endverbatim

#include <vcl_vector.h>
#include <vcl_string.h>
#include "dbdet_edgemap_sptr.h"
#include <vsol/vsol_spatial_object_2d_sptr.h>

//: Loads an ascii file describing edgels.
//\author Amir Tamrakar
bool dbdet_load_edg(vcl_string input_file, bool bSubPixel, bool blines, double scale,
                    vcl_vector< vsol_spatial_object_2d_sptr > &edgels);

bool dbdet_load_edg(vcl_string input_file, bool bSubPixel, double scale, dbdet_edgemap_sptr &edge_map);

//: Save an edgemap as a .edg file
//\author Amir Tamrakar
bool dbdet_save_edg(vcl_string filename, dbdet_edgemap_sptr edgemap);

#endif // dbdet_load_edg_h

