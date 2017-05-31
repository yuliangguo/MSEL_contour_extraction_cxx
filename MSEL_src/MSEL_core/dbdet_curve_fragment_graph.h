// This is dbdet_curve_fragment_graph.h
#ifndef dbdet_curve_fragment_graph_h
#define dbdet_curve_fragment_graph_h
//:
//\file
//\brief Curve Fragment graph structure
//\author Amir Tamrakar
//\date 04/20/07
//
//\verbatim
//  Modifications
//
//    Amir Tamrakar  The graph structure used to be represented as a node adjacency list
//                   where the nodes are stored in a vector indexed by its id. This was simply
//                   mimicking the edge link graph setup. However, this is a waste for representing
//                   a curve fragment graph where the nodes are only the endpoint and junction edgels.
//                   So I've switched to maps instead.
//\endverbatim

#include <vcl_vector.h>
#include <vcl_list.h>
#include <vcl_set.h>
#include "dbdet_edgel.h"
#include "dbdet_edgemap.h"
#include "dbdet_CFTG.h"

//: This class represents the curve fragment graph formed from the edgels
//  The links are curve fragments represented by edgel chains.
class dbdet_curve_fragment_graph
{
public:
  vcl_vector<dbdet_edgel_chain_list> cFrags; ///< child curve fragments
  vcl_vector<dbdet_edgel_chain_list> pFrags; ///< parent curve fragments

  dbdet_edgel_chain_list frags; ///< redundant single list of all fragments
  
  //this is a hack (need to move this out of here into the storage class)
  dbdet_CFTG CFTG; ///< The Curve Fragment Topology Graph (CFTG) 

  vcl_set<int> participate_edge_id;  // include edges in unambiguous fragments and edges participate in hypothesis tree;

  //: constructor
  dbdet_curve_fragment_graph(int size=0): cFrags(size), pFrags(size){}

  //: destructor
  ~dbdet_curve_fragment_graph()
  {
    clear(); //delete everything upon exit
    CFTG.clear();
  }

  //Access functions
  unsigned size() { return cFrags.size(); }//should be the same as edgels.size()

  //: resize the graph
  void resize(unsigned size)
  { 
    if (size!=cFrags.size()){
      clear();
      CFTG.clear();
    }

    cFrags.resize(size);
    pFrags.resize(size);

    CFTG.resize(size);
  }

  //: clear the graph
  void clear()
  {
    //delete all the curve fragments
    dbdet_edgel_chain_list_iter f_it = frags.begin();
    for (; f_it != frags.end(); f_it++)
      delete (*f_it);

    frags.clear();
    cFrags.clear();
    pFrags.clear();

    CFTG.clear();
  }

  //: add a curve fragment to the graph
  void insert_fragment(dbdet_edgel_chain* chain)
  {
    dbdet_edgel* e1 = chain->edgels.front();
    dbdet_edgel* e2 = chain->edgels.back();

    cFrags[e1->id].push_back(chain);
    pFrags[e2->id].push_back(chain);

    frags.push_back(chain);
  }

  //: remove a curve fragment
  void remove_fragment(dbdet_edgel_chain* chain)
  {
    dbdet_edgel* e1 = chain->edgels.front();
    dbdet_edgel* e2 = chain->edgels.back();

    pFrags[e2->id].remove(chain);
    cFrags[e1->id].remove(chain);

    frags.remove(chain);

    delete chain;
  }

  //: just extract a curve fragment from the graph do not delete
  void extract_fragment(dbdet_edgel_chain* chain)
  {
    dbdet_edgel* e1 = chain->edgels.front();
    dbdet_edgel* e2 = chain->edgels.back();

    pFrags[e2->id].remove(chain);
    cFrags[e1->id].remove(chain);

    frags.remove(chain);
  }

  friend class dbdet_edge_map;
};

#endif // dbdet_curve_fragment_graph_h
