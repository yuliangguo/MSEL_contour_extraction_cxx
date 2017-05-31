//This is brcv/seg/dbdet/pro/dbdet_load_edg_process.cxx

#include <vcl_iostream.h>
#include <vcl_cassert.h>
#include <vcl_fstream.h>
#include <vcl_cmath.h>
#include <vcl_algorithm.h>
#include <vul/vul_file.h>
#include <vul/vul_file_iterator.h>

#include "dbdet_load_edg_process.h"

#include <vsol/vsol_point_2d.h>
#include <vsol/vsol_line_2d.h>
#include <vsol/vsol_line_2d_sptr.h>
#include <vgl/vgl_point_2d.h>
#include <vgl/vgl_vector_2d.h>

//#include "dbdet_config.h"
#include "dbdet_edgemap_storage.h"
#include "dbdet_edgemap_storage_sptr.h"

#include "dbdet_edgemap.h"
#include "dbdet_edgemap_sptr.h"
#include "dbdet_load_edg.h"

dbdet_load_edg_process::dbdet_load_edg_process() : bpro1_process(), num_frames_(0)
{
  if( 
#ifdef HAS_BOOST
      !parameters()->add( "Input file <filename.edg..> (gzipped .gz supported)" , "-edginput" , bpro1_filepath("","*.edg*") ) ||
#else
      !parameters()->add( "Input file <filename.edg..> (gzipped .gz NOT supported)" , "-edginput" , bpro1_filepath("","*.edg") ) ||
#endif
      !parameters()->add( "SubPixel ?" ,               "-bP_SP" ,    true ) ||
      !parameters()->add( "Load as vsol" ,             "-bvsol" ,    false )  ||
      !parameters()->add( "Load as Lines" ,            "-blines" ,   true )  ||
      !parameters()->add( "Scale" ,                    "-scale" ,   1.0 ) ||
      !parameters()->add( "Order by filename"        , "-orderbf", false)
    )
  {
    vcl_cerr << "ERROR: Adding parameters in " __FILE__ << vcl_endl;
  }
}

vcl_string dbdet_load_edg_process::name() 
{
  return "Load .EDG File";
}

vcl_vector< vcl_string > dbdet_load_edg_process::get_input_type() 
{
  vcl_vector< vcl_string > to_return;
  return to_return;
}

vcl_vector< vcl_string > dbdet_load_edg_process::get_output_type() 
{
  vcl_vector< vcl_string > to_return;
  bool bvsol;
  parameters()->get_value( "-bvsol" , bvsol );

  if (bvsol)
    to_return.push_back( "vsol2D" );
  else
    to_return.push_back( "edge_map");

  return to_return;
}

//: Clone the process
bpro1_process*
dbdet_load_edg_process::clone() const
{
  return new dbdet_load_edg_process(*this);
}


bool dbdet_load_edg_process::execute()
{
  bpro1_filepath input;
  parameters()->get_value( "-edginput" , input);
  vcl_string input_file_path = input.path;

  int num_of_files = 0;

  output_data_.clear();

  // make sure that input_file_path is sane
  if (input_file_path == "") { return false; }

  //vcl_cout << vul_file::dirname(input_file_path);

  // test if fname is a directory
  if (vul_file::is_directory(input_file_path))
  {
    vcl_vector<vcl_string> input_files;
    for(vul_file_iterator fn = input_file_path + "/*.edg*"; fn; ++fn)
      input_files.push_back(fn());

    // Sort - because the file iterator uses readdir() it does not
    //        iterate over files in alphanumeric order
    bool orderbf;
    parameters()->get_value( "-orderbf" , orderbf );
    if(orderbf) 
      vcl_sort(input_files.begin(), input_files.end());

    while (!input_files.empty())
    {
      vcl_string input_file = input_files.back();
      loadEDG(input_file);
      num_of_files++;
      input_files.pop_back();
    }

    //this is the number of frames to be outputted
    num_frames_ = num_of_files;
  }
  else {
    vcl_string input_file = input_file_path;
    bool successful = loadEDG(input_file);
    num_frames_ = 1;

    if (!successful)
      return false;
  }

  //reverse the order of the objects so that they come out in the right order
  // vcl_reverse(output_data_.begin(),output_data_.end());

  return true;
}

bool dbdet_load_edg_process::loadEDG(vcl_string input_file)
{
  bool bSubPixel=false, blines=false, bvsol=false;
  double scale;
  
  parameters()->get_value( "-bP_SP" , bSubPixel );
  parameters()->get_value( "-bvsol" , bvsol );
  parameters()->get_value( "-blines" , blines );
  parameters()->get_value( "-scale" , scale );

  if (bvsol){
    // edgels (vsol)
    vcl_vector< vsol_spatial_object_2d_sptr > edgels;

    bool retval = dbdet_load_edg(input_file, bSubPixel, blines, scale, edgels);
    if (!retval)
      return false;

    vcl_cout << "N edgels: " << edgels.size() << vcl_endl;
    // create the output storage class
    vidpro1_vsol2D_storage_sptr output_vsol = vidpro1_vsol2D_storage_new();
    output_vsol->add_objects(edgels, input_file);
    output_data_.push_back(vcl_vector< bpro1_storage_sptr > (1,output_vsol));
  }
  else {
    // edge_map 
    dbdet_edgemap_sptr edge_map;

    bool retval = dbdet_load_edg(input_file, bSubPixel, scale, edge_map);
    if (!retval)
      return false;

    vcl_cout << "N edgels: " << edge_map->num_edgels() << vcl_endl;
    // create the output storage class
    dbdet_edgemap_storage_sptr output_edgemap = dbdet_edgemap_storage_new();
    output_edgemap->set_edgemap(edge_map);
    output_data_.push_back(vcl_vector< bpro1_storage_sptr > (1,output_edgemap));
  }
 
  vcl_cout << "Loaded: " << input_file.c_str() << ".\n";

  return true;
}

