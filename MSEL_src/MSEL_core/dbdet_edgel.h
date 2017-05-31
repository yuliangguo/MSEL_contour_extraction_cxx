// This is dbdet_edgel.h
#ifndef dbdet_edgel_h
#define dbdet_edgel_h
//:
//\file
//\brief Edgel class for the Symbolic edge linking algorithm
//\author Amir Tamrakar
//\date 09/05/06
//
//\verbatim
//  Modifications
//  Amir Tamrakar 09/05/06          Moved it from dbdet_se1.h to a new file
//  Amir Tamrakar 11/15/06          Removed the link graph structures from this class to make it lighter
//
//  Ozge Can Ozcanli Jan 12, 2007   Added copy constructor
//  Ricardo Fabbri Aug 03, 2009     Fixed copy constructor to clone appearance members, 
//                                  to avoid delete-ing unitialized memory.
//  Ricardo Fabbri Aug 12, 2009     Added copy assignment, avoiding segfaults.
//
//\endverbatim

#include <vbl/vbl_ref_count.h>
#include <vcl_vector.h>
#include <vcl_deque.h>
#include <vcl_list.h>

#include <vgl/vgl_point_2d.h>

#include "dbdet_appearance.h"

//forward class definitions
class dbdet_edgel;
class dbdet_curvelet;

//useful type definitions
typedef vcl_deque<dbdet_edgel* > dbdet_edgel_list;
typedef vcl_deque<dbdet_edgel* >::iterator dbdet_edgel_list_iter;
typedef vcl_deque<dbdet_edgel* >::const_iterator dbdet_edgel_list_const_iter;
typedef vcl_deque<dbdet_edgel* >::reverse_iterator dbdet_edgel_list_reverse_iter;
typedef vcl_deque<dbdet_edgel* >::const_reverse_iterator dbdet_edgel_list_const_reverse_iter;

typedef vcl_list<dbdet_curvelet* > curvelet_list;
typedef vcl_list<dbdet_curvelet* >::iterator curvelet_list_iter;
typedef vcl_list<dbdet_curvelet* >::const_iterator curvelet_list_const_iter;

//: edgel class: contains pt, tangent and collection of all the groupings around it
class dbdet_edgel
{
public:
  int id;                  ///< unique id

  vgl_point_2d<double> pt; ///< the location of the edgel
  double tangent;          ///< the orientation of the edgel in radians

  //: the strength of the edgel (typically gradient magnitude)
  double strength;         

  //: second derivative of the image in the direction of the gradient
  double deriv;            

  //: Harris Corner like measure to factor in the error in edge localization at
  // corners and junctions
  double uncertainty;      
  
  vgl_point_2d<int> gpt;   ///< hack to store the grid coordinates

  //: appearance information stored with the edgels
  dbdet_appearance* left_app;  ///< appearance on the left side of the edgel
  dbdet_appearance* right_app; ///< appearance on the right side of the edgel

  //: default constructor
  dbdet_edgel() : 
    id((unsigned)-1), 
    pt(vgl_point_2d<double>(0,0)), 
    tangent(0.0), 
    strength(0.0), 
    deriv(0.0), 
    uncertainty(0.0), 
    left_app(new dbdet_intensity(0.0)), 
    right_app(new dbdet_intensity(0.0)){}
  
  //: constructor
  //
  dbdet_edgel(vgl_point_2d<double> new_pt, double tan, 
      double conf=0.0, double der=0.0, double uncer=0.0,
      dbdet_appearance* lapp=new dbdet_intensity(0.0), 
      dbdet_appearance* rapp=new dbdet_intensity(0.0));

  //: copy constructor
  dbdet_edgel(const dbdet_edgel& other);

  //: Assignment operator
  dbdet_edgel & operator=(const dbdet_edgel &rhs);

  //: destructor
  ~dbdet_edgel();
};

//: A class to hold a chain of edgels that defines the image curve
class dbdet_edgel_chain
{
public:
  dbdet_edgel_list edgels;
  bool temp; //temp flag for the CFTG 

  //: constructor
  dbdet_edgel_chain(): edgels(0), temp(false){}
  ~dbdet_edgel_chain(){}

  //: copy constructor
  dbdet_edgel_chain(const dbdet_edgel_chain& chain): 
    edgels(chain.edgels.size()), temp(false)
  {
    for (unsigned i=0; i<chain.edgels.size(); i++)
      edgels[i] = chain.edgels[i];
  }

  //merge two edge chains together (expect some overlap)
  void merge(dbdet_edgel_chain* n_chain)
  {
    //there are going to be some duplications 
    //need to remove those elements
    bool found_last=false;
    dbdet_edgel_list_iter lit = n_chain->edgels.begin();
    for(;lit!=n_chain->edgels.end();lit++){
      if ((*lit)==edgels.back())
        found_last=true;
      
      if (found_last)
        edgels.push_back(*lit);
    }
  }

  void append(vcl_vector<dbdet_edgel*>& n_chain)
  {
    for (unsigned i=0; i<n_chain.size(); i++)
      edgels.push_back(n_chain[i]);
  }

  void append(dbdet_edgel_list& n_chain)
  { 
    int sizeofchain=n_chain.size();
    for (unsigned i=1; i<sizeofchain; i++)
      edgels.push_back(n_chain[i]);
  }

  void push_back(dbdet_edgel* edgel){ edgels.push_back(edgel); }
  void push_front(dbdet_edgel* edgel){ edgels.push_front(edgel); }
};

typedef vcl_list<dbdet_edgel_chain*> dbdet_edgel_chain_list;
typedef vcl_list<dbdet_edgel_chain*>::iterator dbdet_edgel_chain_list_iter;
typedef vcl_list<dbdet_edgel_chain*>::const_iterator dbdet_edgel_chain_list_const_iter;

#endif // dbdet_edgel_h
