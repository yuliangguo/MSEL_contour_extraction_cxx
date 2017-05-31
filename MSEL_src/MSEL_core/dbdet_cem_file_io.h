// This is brcv/seg/dbdet/algo/dbdet_cem_file_io.h
#ifndef dbdet_cem_file_io_h
#define dbdet_cem_file_io_h
//:
//\file
//\brief Utility functions for contour file load/save
//\author Amir Tamrakar
//\date November 15, 2007
//
//\verbatim
//Modifications
// firat: cem_v1 can now be loaded with correct width and height instead of using the default value 1000.
// firat: converting degrees to radians is now optional while loading cem_v1 files.
//
// \endverbatim

#include <vcl_fstream.h>
#include <vcl_vector.h>
#include <vcl_string.h>
#include "dbdet_edgemap_sptr.h"
#include "dbdet_curve_fragment_graph.h"

//: Save the contour fragment graph as a .cem file
bool dbdet_save_cem(vcl_string filename, dbdet_edgemap_sptr EM, dbdet_curve_fragment_graph& CFG);

//: Loads an ascii file containing a graph of edgel chains (the contour fragment graph)
dbdet_edgemap_sptr dbdet_load_cem(vcl_string filename, dbdet_curve_fragment_graph& CFG);

//: load cem file version 2
dbdet_edgemap_sptr dbdet_load_cem_v2(vcl_ifstream &infp, dbdet_curve_fragment_graph& CFG);

//: load the older version (aka version 1) of the cem file (for backward compatibility)
dbdet_edgemap_sptr dbdet_load_cem_v1(vcl_ifstream &infp, dbdet_curve_fragment_graph& CFG, int width = 3000, int height = 3000, bool convert_degrees_to_radians = true);

#endif // dbdet_cem_file_io_h

