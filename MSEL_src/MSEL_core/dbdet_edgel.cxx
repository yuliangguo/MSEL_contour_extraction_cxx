#include "dbdet_edgel.h"
#include "dbdet_sel_utils.h"

#include <vcl_iostream.h>
#include <vcl_fstream.h>
#include <vcl_cassert.h>
#include <vcl_deque.h>
#include <vcl_algorithm.h>

//------------------------------------------------------------------------------
// dbdet_edgel methods
//------------------------------------------------------------------------------

//: constructor
dbdet_edgel::dbdet_edgel(
    vgl_point_2d<double> new_pt, double tan, double conf, double der, double uncer, 
    dbdet_appearance* lapp, dbdet_appearance* rapp) 
: 
  id(-1), 
  pt(new_pt), 
  tangent(dbdet_angle0To2Pi(tan)), 
  strength(conf), 
  deriv(der), 
  uncertainty(uncer), 
  left_app(lapp), 
  right_app(rapp)
{
}


//: copy constructor
dbdet_edgel::dbdet_edgel(const dbdet_edgel& other)
  :
  id(other.id),
  pt(other.pt),
  tangent(other.tangent),
  strength(other.strength),
  deriv(other.deriv),
  gpt(other.gpt)
{
  left_app = other.left_app->clone();
  right_app = other.right_app->clone();
}

dbdet_edgel & 
dbdet_edgel::
operator=(const dbdet_edgel &rhs)
{
  id = rhs.id;

  pt = rhs.pt;
  tangent = rhs.tangent;
  strength = rhs.strength;
  deriv = rhs.deriv;
  
  gpt = rhs.gpt;

  //: Copy pointers, taking care for when rhs is *this
  dbdet_appearance* left_app_orig = left_app;
  dbdet_appearance* right_app_orig = right_app;

  left_app = rhs.left_app->clone();
  right_app = rhs.right_app->clone();

  delete left_app_orig;
  delete right_app_orig;

  return *this;
}

//: destructor
dbdet_edgel::~dbdet_edgel()
{
  delete left_app;
  delete right_app;
}
