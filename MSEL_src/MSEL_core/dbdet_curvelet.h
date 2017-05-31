// This is dbdet_curvelet.h
#ifndef dbdet_curvelet_h
#define dbdet_curvelet_h
//:
//\file
//\brief A class to represent an edgel grouping with an associated curve model (curvelet now called curve bundle)
//\author Amir Tamrakar
//\date 09/05/06
//
//\verbatim
//  Modifications
//  Amir Tamrakar 09/05/06          Moved it from dbdet_se1.h to a new file
//
//  Ozge Can Ozcanli Jan 12, 2007   Added copy constructor
//
//\endverbatim

#include <vbl/vbl_ref_count.h>
#include <vcl_vector.h>
#include <vcl_list.h>
#include <vcl_set.h>
#include <vcl_map.h>
#include <vcl_utility.h>

#include <vnl/vnl_math.h>

#include <vgl/vgl_point_2d.h>
#include <vgl/vgl_polygon.h>
#include <vgl/vgl_distance.h>
#include <vgl/vgl_clip.h>
#include <vgl/vgl_area.h>

#include "bgld_eulerspiral.h"
#include "dbdet_edgel.h"
#include "dbdet_curve_model.h"

//forward declaration
class dbdet_curvelet;

typedef vcl_list<dbdet_curvelet*> cvlet_list;
typedef vcl_list<dbdet_curvelet*>::iterator cvlet_list_iter;

//: The curvelet class: stores the ordered list of edgels defining the curvelet
//  and also the curve model defined by the grouping
//  It also stores a list of higher order curvelets that it forms.
class dbdet_curvelet
{
public:
  dbdet_edgel* ref_edgel;                        ///< ref edgel (the edgel to which it is anchored)

  vcl_vector<dbdet_edgel*> edgel_chain;          ///< the ordered list of edgels 
  dbdet_curve_model* curve_model;                ///< associated curve model

  bool forward;                                  ///< is this a forward cvlet or a reverse cvlet
  double length;                                 ///< length of the curvelet
  double quality;                                ///< the quality of this grouping (determined by various means)
  
  bool used;                                     ///< to keep track of whether this curvelet was used in linking

  //curvelet_list component_curvelets;             ///< list of curvelets forming this curvelet
  //curvelet_list larger_curvelets;                ///< list of larger curvelets formed by this curvelet

  //: default constructor
  dbdet_curvelet(): ref_edgel(0), edgel_chain(0), curve_model(0), forward(true), length(0.0), quality(0.0), used(false){}

  //: constructor 1
  dbdet_curvelet(dbdet_edgel* e) : ref_edgel(e), edgel_chain(0), curve_model(0), forward(true), length(0.0), quality(0.0), used(false){}

  //: constructor 2
  dbdet_curvelet(dbdet_edgel* e, dbdet_curve_model* cm, bool dir=true) : 
    ref_edgel(e), edgel_chain(0), curve_model(cm), forward(dir), length(0.0), quality(0.0), used(false){}

  //: constructor 3
  dbdet_curvelet(dbdet_edgel* e, dbdet_curve_model* cm, vcl_vector<dbdet_edgel*> &echain, bool dir=true) : 
    ref_edgel(e), edgel_chain(echain), curve_model(cm), forward(dir), length(0.0), quality(0.0), used(false){}

  //: constructor 4
  dbdet_curvelet(dbdet_edgel* e, dbdet_curve_model* cm, vcl_deque<dbdet_edgel*> &echain, bool dir=true) : 
    ref_edgel(e), curve_model(cm), forward(dir), length(0.0), quality(0.0), used(false)
  {
    edgel_chain.insert(edgel_chain.end(), echain.begin(), echain.end());
  }

  //: copy constructor
  dbdet_curvelet(const dbdet_curvelet& other);

  //: copy constructor with provisions for copying a different curve bundle
  dbdet_curvelet(const dbdet_curvelet& other, dbdet_curve_model* cm);

  //: destructor
  ~dbdet_curvelet();

  //: return the order of this grouping
  unsigned order() const { return edgel_chain.size(); }

  //: add an edgel 
  void push_back(dbdet_edgel* e) { edgel_chain.push_back(e); }

  //: add the associated curve model for this grouping
  void set_curve_model(dbdet_curve_model* cm) { curve_model = cm; }

  //: replace the curve model associated with this grouping with a new one
  void replace_curve_model(dbdet_curve_model* cm) 
  { 
    if (!cm){ 
      vcl_cout << "attempting to replace CB with NULL!" << vcl_endl; 
      return; 
    }

    //assert(cm); 
    if (curve_model) 
      delete curve_model; 
    
    curve_model = cm; 
  }

  //: returns the forward(child) edgel_chain 
  vcl_list<dbdet_edgel*> child_chain();

  //: returns the backward(child) edgel_chain 
  vcl_list<dbdet_edgel*> parent_chain();

  ////: add a larger grouping
  //void add_larger_curvelet(dbdet_curvelet* cv) { larger_curvelets.push_back(cv); }

  ////: does a larger grouping exist already?
  //dbdet_curvelet* exists(dbdet_curvelet* cv);

  //: compute properties of this curvelet once formed
  void compute_properties(double R, double token_len);

  //: print info to file
  void print(vcl_ostream&);

  //: read info from file
  void read(vcl_istream&);

};

#endif // dbdet_curvelet_h
