#include "dbdet_sel_base.h"

#include <vcl_iostream.h>
#include <vcl_fstream.h>
#include <vcl_cassert.h>
#include <vcl_deque.h>
#include <vcl_map.h>
#include <vcl_set.h>
#include <vcl_algorithm.h>
#include <vcl_queue.h>



//********************************************************************//
// Functions for constructing hypothesis trees and
// tracing between contour end points.
//********************************************************************//


#define SM_TH 0.6

#define Theta_1 0.6
#define Theta_2 0.0
#define Strength_Diff1 2.0
#define Strength_Diff2 4.0

dbdet_EHT* dbdet_sel_base::construct_hyp_tree(dbdet_edgel* edge)
{
  if (edge_link_graph_.cLinks.size()==0){
    vcl_cout << "No Link Graph !" <<vcl_endl;
    return 0;
  }


 //construct 2 HTs: one in the forward direction and one in the reverse direction ????  by yuliang no forword
//Modify: for each node, consider all the child links and parent links(exact the one linking parent node), as its child node
  vcl_queue<dbdet_EHT_node*> BFS_queue;

  //forward HT
  dbdet_EHT* HTF = new dbdet_EHT();

  dbdet_EHT_node* root1 = new dbdet_EHT_node(edge);
  HTF->root = root1;
  BFS_queue.push(root1);

  int depth = 0; // comment by Yuliang, this is not the depth of the tree, but number of nodes actually

  //How far do we wanna go (if we don't hit a node)?
  while (!BFS_queue.empty() && vcl_log10(double(depth))<3)
  {
    dbdet_EHT_node* cur_node = BFS_queue.front();
    BFS_queue.pop();

    //are we at a CFG node? if we are we don't need to go any further
    if (cur_node!= root1 &&
        (curve_frag_graph_.pFrags[cur_node->e->id].size()>0 ||
         curve_frag_graph_.cFrags[cur_node->e->id].size()>0))
      continue;

    //also if we hit an edgel that is already linked no need to go further (this might ensure planarity)
    if (edge_link_graph_.linked[cur_node->e->id])
      continue;

    //propagate this node
    dbdet_link_list_iter lit = edge_link_graph_.cLinks[cur_node->e->id].begin();
    for (; lit != edge_link_graph_.cLinks[cur_node->e->id].end(); lit++)
    {
      if (edge_link_graph_.linked[(*lit)->ce->id]) //don't go tracing in linked contours
        continue;
  
      if (cur_node->parent) {
        //make a simple consistency check
        double dx1 = cur_node->e->pt.x() - cur_node->parent->e->pt.x();
        double dy1 = cur_node->e->pt.y() - cur_node->parent->e->pt.y();
        double dx2 = (*lit)->ce->pt.x() - cur_node->e->pt.x();
        double dy2 = (*lit)->ce->pt.y() - cur_node->e->pt.y();


        if (((dx1*dx2 + dy1*dy2)/vcl_sqrt(dx1*dx1+dy1*dy1)/vcl_sqrt(dx2*dx2+dy2*dy2))<SM_TH) //not consistent, but with lower TH to keep more Hypothesis by Yuliang ////////// Cosine Formula
          continue;
      }

      //else extend the tree to this edgel
      dbdet_EHT_node* new_node = new dbdet_EHT_node((*lit)->ce);

      cur_node->add_child(new_node);
      BFS_queue.push(new_node);
      depth++;
    }
    // by Yuliang
    // explore in pLinks
    lit = edge_link_graph_.pLinks[cur_node->e->id].begin();
    for (; lit != edge_link_graph_.pLinks[cur_node->e->id].end(); lit++)
    {
      if (edge_link_graph_.linked[(*lit)->pe->id]) //don't go tracing in linked contours
        continue;
      if (cur_node->parent) {
          if((*lit)->pe->id == cur_node->parent->e->id)// if this parent link the same from parent node, don't trace
        continue;
        //make a simple consistency check
        double dx1 = cur_node->e->pt.x() - cur_node->parent->e->pt.x();
        double dy1 = cur_node->e->pt.y() - cur_node->parent->e->pt.y();
        double dx2 = (*lit)->pe->pt.x() - cur_node->e->pt.x();
        double dy2 = (*lit)->pe->pt.y() - cur_node->e->pt.y();

        if (((dx1*dx2 + dy1*dy2)/vcl_sqrt(dx1*dx1+dy1*dy1)/vcl_sqrt(dx2*dx2+dy2*dy2))<SM_TH) //not consistent, but with lower TH to keep more Hypothesis by Yuliang
          continue;
      }

      //else extend the tree to this edgel
      dbdet_EHT_node* new_node = new dbdet_EHT_node((*lit)->pe);

      cur_node->add_child(new_node);
      BFS_queue.push(new_node);
      depth++;
    }
  }

  //empty the bfs queue
  while (!BFS_queue.empty())
    BFS_queue.pop();

  return HTF;
}

//: construct all possible EHTs from the terminal nodes and find legal contour paths
void dbdet_sel_base::construct_all_path_from_EHTs()
{
  // modify by Yuliang, use a local filter for contours size 3 and 4 to prune some, and merge some continuous regular contours first
  regular_contour_filter();
  //go over the contour fragment graph and form an EHT from every terminal node
  //validate each of the paths in the EHT

  vcl_vector<dbdet_edgel_chain*> new_frags;

  //going over the edgemap instead so that an EHT only starts once from a node when there are two
  //contour fragments terminating there
  for (unsigned i=0; i<edgemap_->edgels.size(); i++)
  {
    dbdet_edgel* eA = edgemap_->edgels[i];

    if (curve_frag_graph_.pFrags[i].size()==0 &&
        curve_frag_graph_.cFrags[i].size()==0)
      continue; //no terminal nodes here

    //1) Terminal node found, construct an EHT from here
    dbdet_EHT* EHT1 = construct_hyp_tree(eA);

    //2) traverse each path and determine its validity
    if (EHT1)
    {
      //traverse the EHT and test all the paths
      dbdet_EHT::path_iterator pit = EHT1->path_begin();
      for (; pit != EHT1->path_end(); pit++){
        vcl_vector<dbdet_edgel*>& edgel_chain = pit.get_cur_path();

    dbdet_edgel* le = edgel_chain.back();
 
        if (curve_frag_graph_.pFrags[le->id].size()==0 && curve_frag_graph_.cFrags[le->id].size()==0)
        {
          //not a valid termination node
          //delete the node associated  with this path ( it will delete the entire path, by definition)
          EHT1->delete_subtree(pit);
          continue;
        }

        //test this path to see if it is valid
        if (!is_EHT_path_legal(edgel_chain)){
          EHT1->delete_subtree(pit);
          continue;
        }

        //copy this chain
        dbdet_edgel_chain* new_chain = new dbdet_edgel_chain();
        new_chain->append(edgel_chain);
        new_chain->temp = true; //make sure that the new frags are initialized with temp flags

	curve_frag_graph_.CFTG.insert_fragment(new_chain);
        //now that its added delete the apth
        EHT1->delete_subtree(pit);
      }

      //finally delete the EHT
      delete EHT1;
    }
  }
  vcl_cout<<"Finish constructing all hypothesis trees"<<vcl_endl;
  ////Now add all the new curve fragments into the CFG (as tentative fragments)
  //for (unsigned i=0; i<new_frags.size(); i++)
  //  curve_frag_graph_.insert_fragment(new_frags[i]);
}

//: perform a geometric consistency check to determine whether a given temp path is valid
bool dbdet_sel_base::is_EHT_path_legal(vcl_vector<dbdet_edgel*>& edgel_chain)
{
  //what makes a path legal?

  // Modified by Yuliang, before (a), (b), if path only have < 3 edgels, prune out direct link over 3 pixels
  if(edgel_chain.size()<3)
  {
      dbdet_edgel* eS = edgel_chain.front();
      dbdet_edgel* eE = edgel_chain.back();
      double dx1 = eE->pt.x() - eS->pt.x();
      double dy1 = eE->pt.y() - eS->pt.y();
      double dist= vcl_sqrt(dx1*dx1+dy1*dy1);
      if(dist>5)
        return false;
   }

  // by Yuliang, construct two lists of end nodes from other Curve Fragments end poits linking the end points of the path
  vcl_vector<dbdet_edgel*> S_link_end_nodes, E_link_end_nodes;

  // (a) if a c1 polyarc bundle can form within it
  // (b) if it is c1 compatible with the end points

  //For now, just check for...

  //1) continuity consistency at the end points
  //1a) at the start point
  dbdet_edgel* eS = edgel_chain.front();
  dbdet_edgel* e2 = edgel_chain[1];
  double dx1 = e2->pt.x() - eS->pt.x();
  double dy1 = e2->pt.y() - eS->pt.y();
  bool cons = false;


  //is this at a start or an end of an unambiguous chain?
  dbdet_edgel_chain_list_iter pcit = curve_frag_graph_.pFrags[eS->id].begin();
  for ( ; pcit != curve_frag_graph_.pFrags[eS->id].end(); pcit++)
  {
    dbdet_edgel* pe = (*pcit)->edgels[(*pcit)->edgels.size()-2];
    S_link_end_nodes.push_back((*pcit)->edgels.front());
 
    //make a simple consistency check for child Frags
    double dx2 = eS->pt.x() - pe->pt.x();
    double dy2 = eS->pt.y() - pe->pt.y();

    cons = cons || ((dx1*dx2 + dy1*dy2)/vcl_sqrt(dx1*dx1+dy1*dy1)/vcl_sqrt(dx2*dx2+dy2*dy2))>0;
  }
  dbdet_edgel_chain_list_iter ccit = curve_frag_graph_.cFrags[eS->id].begin();
  for ( ; ccit != curve_frag_graph_.cFrags[eS->id].end(); ccit++)
  {
    dbdet_edgel* ce = (*ccit)->edgels[1];
    S_link_end_nodes.push_back((*ccit)->edgels.back());

    //make a simple consistency check for parent Frags
    double dx2 = eS->pt.x() - ce->pt.x();
    double dy2 = eS->pt.y() - ce->pt.y();

    cons = cons || ((dx1*dx2 + dy1*dy2)/vcl_sqrt(dx1*dx1+dy1*dy1)/vcl_sqrt(dx2*dx2+dy2*dy2))>0;
  }
  if (!cons) return false; //no good at the start point

  //1b) at the end point
  dbdet_edgel* eE = edgel_chain.back();
  e2 = edgel_chain[edgel_chain.size()-2];
  dx1 = eE->pt.x() - e2->pt.x();
  dy1 = eE->pt.y() - e2->pt.y();
  cons = false;
 
  //is this at a start or an end of an unambiguous chain?
  pcit = curve_frag_graph_.pFrags[eE->id].begin();
  for ( ; pcit != curve_frag_graph_.pFrags[eE->id].end(); pcit++)
  {
    dbdet_edgel* pe = (*pcit)->edgels[(*pcit)->edgels.size()-2];
    E_link_end_nodes.push_back((*pcit)->edgels.front());

    //make a simple consistency check
    double dx2 = pe->pt.x() - eE->pt.x();
    double dy2 = pe->pt.y() - eE->pt.y();

    cons = cons || ((dx1*dx2 + dy1*dy2)/vcl_sqrt(dx1*dx1+dy1*dy1)/vcl_sqrt(dx2*dx2+dy2*dy2))>0;
  }
  ccit = curve_frag_graph_.cFrags[eE->id].begin();
  for ( ; ccit != curve_frag_graph_.cFrags[eE->id].end(); ccit++)
  {
    dbdet_edgel* ce = (*ccit)->edgels[1];
    E_link_end_nodes.push_back((*ccit)->edgels.back());

    //make a simple consistency check
    double dx2 = ce->pt.x() - eE->pt.x();
    double dy2 = ce->pt.y() - eE->pt.y();

    cons = cons || ((dx1*dx2 + dy1*dy2)/vcl_sqrt(dx1*dx1+dy1*dy1)/vcl_sqrt(dx2*dx2+dy2*dy2))>0;
  }
  if (!cons) return false; //no good at the end point

  //2) By Yuliang, use the two lists, check if it is a path between which unambiguous contours linked
  for(int i =0; i<S_link_end_nodes.size(); i++)
  {
      dbdet_edgel* S_l_e = S_link_end_nodes[i];
    if(S_l_e == eE)
        return false;
      for(int j=0; j<E_link_end_nodes.size(); j++)
        {
        dbdet_edgel* E_l_e = E_link_end_nodes[j];
        if(E_l_e == eS)
            return false;
        // the case two connected contours filling the path
        if(E_l_e== S_l_e)
            return false;
        }
  }

  // comment by Yuliang, in most cases, it is no need, because of short paths.
  //fit_polyarc_to_chain(&edgel_chain);
  return true;
}

//: New Quality Metric by Naman Kumar :: compute a path metric based on the Gap, Orientation, Strength and Size of the chain
double dbdet_sel_base::compute_path_metric2(vcl_vector<dbdet_edgel*>& Pchain,
                                           vcl_vector<dbdet_edgel*>& Tchain,
                                           vcl_vector<dbdet_edgel*>& Cchain)
{
  double cost = 0.0;double ds=0;double dt=0;

  //construct an edgel chain out of all three chains
  vcl_vector<dbdet_edgel*> chain;
  if (Pchain.size())
    for (unsigned i=0; i<Pchain.size(); i++) chain.push_back(Pchain[i]);
  if (Tchain.size())
    for (unsigned i=0; i<Tchain.size(); i++) chain.push_back(Tchain[i]);
  if (Cchain.size())
    for (unsigned i=0; i<Cchain.size(); i++) chain.push_back(Cchain[i]);

  //now compute the metric
  dbdet_edgel *eA=0, *eP=0;
  double dsp = 0, thp = 0, total_ds =0.0, a=0.0,s1=0,s2=0,s=0,size=chain.size();
  for (unsigned i=1; i<chain.size(); i++)
  {
    // computing difference in strength
    eA = chain[i];
    eP = chain[i-1];
    s1=(eA)->strength;
    s2=(eP)->strength;
    s=vcl_fabs(s1-s2);
    //computing ds
    ds = vgl_distance(eA->pt, eP->pt);
    if(ds>1.0) a=2.0; else a=1.0;
    total_ds += ds;
    //computing dtheta
    double thc = dbdet_vPointPoint(eP->pt, eA->pt);
    dt = vcl_fabs(thc-thp);
    dt = (dt>vnl_math::pi)? 2*vnl_math::pi-dt : dt;
    cost += vcl_pow((s+dt + a*ds)/size, 2.0); 
    thp = thc;//saving the current vector for the next iteration
    dsp = ds;
  }
  return cost;
}

/*//: New Quality Metric by Yuliang Guo :: compute a path metric based on the Gap, Orientation normorlized by Size of the chain, weights trained by grid search
double dbdet_sel_base::compute_path_metric2(vcl_vector<dbdet_edgel*>& Pchain,
                                           vcl_vector<dbdet_edgel*>& Tchain,
                                           vcl_vector<dbdet_edgel*>& Cchain)
{
  double cost = 0.0;double ds=0;double dt=0;

  //construct an edgel chain out of all three chains
  vcl_vector<dbdet_edgel*> chain;
  if (Pchain.size())
    for (unsigned i=0; i<Pchain.size(); i++) chain.push_back(Pchain[i]);
  if (Tchain.size())
    for (unsigned i=0; i<Tchain.size(); i++) chain.push_back(Tchain[i]);
  if (Cchain.size())
    for (unsigned i=0; i<Cchain.size(); i++) chain.push_back(Cchain[i]);

  //now compute the metric
  dbdet_edgel *eA=0, *eP=0;
  double dsp = 0, thp = 0, total_ds =0.0, w=0.2,s1=0,s2=0,s=0,size=chain.size();
  for (unsigned i=1; i<chain.size(); i++)
  {
    eA = chain[i];
    eP = chain[i-1];
    //s1=(eA)->strength;
    //s2=(eP)->strength;
    //s=vcl_fabs(s1-s2);
    //compute ds
    ds = vgl_distance(eA->pt, eP->pt);
    //if(ds>1.0) a=2.0; else a=1.0;
    total_ds += ds;
    //compute dtheta
    double thc = dbdet_vPointPoint(eP->pt, eA->pt);
    dt = vcl_fabs(thc-thp);
    dt = (dt>vnl_math::pi)? 2*vnl_math::pi-dt : dt;
    //cost += vcl_pow((s+dt + a*ds)/size, 2.0); 
    cost += (w*dt + (1-w)*ds)/size; 
    thp = thc;//save the current vector for the next iteration
    dsp = ds;
  }
  return cost;
}*/



bool link_cost_less(dbdet_CFTG_link* link1, dbdet_CFTG_link* link2)
{return link1->cost < link2->cost;}


//: disambiguate the CFG, basically to produce a disjoint set
void dbdet_sel_base::disambiguate_the_CFTG()
{
  //At the moment, I cannot verify that the CFTG is topologically sound (i.e., a planar graph)
  //so within this limit, the goal is to break the links at junctions

  //Alternatively, it is possible to splice the two contours and mark the connection with the others as a junction
  //these others might be pruned off at a postprocessing stage if desired

  //Note: the temp flag on the contours distinguish it from the unambiguous contours

  //A) disambiguate the links first : Only keep the best path
  //   Note: remember to search in both directions

  //go over all the links of the CFTG
  vcl_vector<dbdet_edgel*> dummy_chain;
   
  dbdet_CFTG_link_list_iter l_it = curve_frag_graph_.CFTG.Links.begin();
  for (; l_it != curve_frag_graph_.CFTG.Links.end(); l_it++)
  {
    dbdet_CFTG_link* cur_Link = (*l_it);
 

    //is this an ambiguous link?
    if (cur_Link->cCFs.size()>1)
    {
      //needs disambiguation
      double min_cost = 10000;
      dbdet_edgel_chain* best_chain = 0;
    dbdet_edgel_chain_list_iter f_it = cur_Link->cCFs.begin();
    for(; f_it != cur_Link->cCFs.end(); f_it++)
      {
        dbdet_edgel_chain* edgel_chain = (*f_it);
        vcl_vector<dbdet_edgel*> chain(edgel_chain->edgels.begin(), edgel_chain->edgels.end());

        double path_cost = compute_path_metric2(dummy_chain, chain, dummy_chain);
        if (path_cost < min_cost){
          min_cost = path_cost;
          best_chain = edgel_chain;
        }
 
      }

      //remove all except the best chain
      if (best_chain){
        dbdet_edgel_chain_list_iter f_it = cur_Link->cCFs.begin();
        for(; f_it != cur_Link->cCFs.end(); f_it++)
          if ((*f_it) != best_chain)
            delete (*f_it);

        cur_Link->cCFs.clear();
        cur_Link->cCFs.push_back(best_chain);
        cur_Link->cost = min_cost;
}
    }
    else { //just comptue cost for this path

      dbdet_edgel_chain* edgel_chain = cur_Link->cCFs.front();
      vcl_vector<dbdet_edgel*> chain(edgel_chain->edgels.begin(), edgel_chain->edgels.end());
      cur_Link->cost = compute_path_metric2(dummy_chain, chain, dummy_chain);
    }
  }
  //B) disambiguate between duplicates (A->B vs B->A)
  l_it = curve_frag_graph_.CFTG.Links.begin();
  for (; l_it != curve_frag_graph_.CFTG.Links.end(); l_it++)
  {
    dbdet_CFTG_link* cur_Link = (*l_it);

    if (cur_Link->cCFs.size()==0)
      continue;

    //look for the link from the other direction
    dbdet_CFTG_link_list_iter l_it2 = curve_frag_graph_.CFTG.cLinks[cur_Link->eE->id].begin();
    for (; l_it2 != curve_frag_graph_.CFTG.cLinks[cur_Link->eE->id].end(); l_it2++){
      if ((*l_it2)->eE == cur_Link->eS){
        //duplicate found

        if ((*l_it2)->cCFs.size()==0)
          continue;

        dbdet_edgel_chain* edgel_chain1 = cur_Link->cCFs.front();
        dbdet_edgel_chain* edgel_chain2 = (*l_it2)->cCFs.front();

        vcl_vector<dbdet_edgel*> chain1(edgel_chain1->edgels.begin(), edgel_chain1->edgels.end());
        vcl_vector<dbdet_edgel*> chain2(edgel_chain2->edgels.begin(), edgel_chain2->edgels.end());

        double path_cost1 = compute_path_metric2(dummy_chain, chain1, dummy_chain);
        double path_cost2 = compute_path_metric2(dummy_chain, chain2, dummy_chain);

        if (path_cost1<path_cost2){
          //keep current link and delete the other
          delete edgel_chain2;
          (*l_it2)->cCFs.clear();
          (*l_it2)->cost = 1000;
        }
        else{
          delete edgel_chain1;
          cur_Link->cCFs.clear();
          cur_Link->cost = 1000;
        }
      }
    }
  }

  //C) Gradient descent to prune bifurcations from the CFTG

  //go over list of Links and find any with degree > 1
  //these need to be disambiguated (gradient descent)
  vcl_list<dbdet_CFTG_link*> GD_list;

  //populate the map
  l_it = curve_frag_graph_.CFTG.Links.begin();
  for (; l_it != curve_frag_graph_.CFTG.Links.end(); l_it++)
  {
    //compute degree at each end
    int deg_S = curve_frag_graph_.CFTG.cLinks[(*l_it)->eS->id].size() + curve_frag_graph_.CFTG.pLinks[(*l_it)->eS->id].size();
    int deg_E = curve_frag_graph_.CFTG.cLinks[(*l_it)->eE->id].size() + curve_frag_graph_.CFTG.pLinks[(*l_it)->eE->id].size();

    if (deg_S>1){
      GD_list.push_back(*l_it);
      continue;
    }

    if (deg_E>1){
      GD_list.push_back(*l_it);
      continue;
    }
 
  }

 
  //sort the cost list
  GD_list.sort(link_cost_less);

  //gradient descent
  while (GD_list.size()>0)
  {
    dbdet_CFTG_link* cur_Link = GD_list.front();
    GD_list.pop_front();

    //now remove the other links connected to the end points of this link
    //clinks from eS
    vcl_vector<dbdet_CFTG_link*> links_to_del;
    l_it = curve_frag_graph_.CFTG.cLinks[cur_Link->eS->id].begin();
    for (; l_it != curve_frag_graph_.CFTG.cLinks[cur_Link->eS->id].end(); l_it++){
 
      if ((*l_it) != cur_Link)
    links_to_del.push_back((*l_it));
    }

    //by yuliang, also consider plinks from eS
    l_it = curve_frag_graph_.CFTG.pLinks[cur_Link->eS->id].begin();
    for (; l_it != curve_frag_graph_.CFTG.pLinks[cur_Link->eS->id].end(); l_it++){

      if ((*l_it) != cur_Link)
        links_to_del.push_back((*l_it));
    }

    for (unsigned j=0; j<links_to_del.size(); j++){   
      GD_list.remove(links_to_del[j]);//also remove it from the GD list
      curve_frag_graph_.CFTG.remove_link(links_to_del[j]);
    }
    links_to_del.clear();

    //plinks from eE
    l_it = curve_frag_graph_.CFTG.pLinks[cur_Link->eE->id].begin();
    for (; l_it != curve_frag_graph_.CFTG.pLinks[cur_Link->eE->id].end(); l_it++){

      if ((*l_it) != cur_Link)
        links_to_del.push_back((*l_it));
    }

    //by yuliang, also consider clinks from eE
    l_it = curve_frag_graph_.CFTG.cLinks[cur_Link->eE->id].begin();
    for (; l_it != curve_frag_graph_.CFTG.cLinks[cur_Link->eE->id].end(); l_it++){
 
      if ((*l_it) != cur_Link)
        links_to_del.push_back((*l_it));
    }

    for (unsigned j=0; j<links_to_del.size(); j++){
      GD_list.remove(links_to_del[j]);//also remove it from the GD list
      curve_frag_graph_.CFTG.remove_link(links_to_del[j]);
    }
    links_to_del.clear();
  }

  //D) Add it all back to the CFG (clear the CFTG in the process)
  l_it = curve_frag_graph_.CFTG.Links.begin();
  for (; l_it != curve_frag_graph_.CFTG.Links.end(); l_it++)
  {
    if ((*l_it)->cCFs.size()==1)
      curve_frag_graph_.insert_fragment((*l_it)->cCFs.front());
  }
  curve_frag_graph_.CFTG.clear();
  curve_frag_graph_.CFTG.resize(edgemap_->edgels.size());
   vcl_cout<<"Finish disambiguating the CFTG"<<vcl_endl;
}

// Following part by Yuliang Guo.
static bool is_continue (const dbdet_edgel_chain *c1, const dbdet_edgel_chain *c2); // test genenal/local continuity

static bool is_longer (const dbdet_edgel_chain *c1, const dbdet_edgel_chain *c2){ // whether contour 1 is longer
    if (c1->edgels.size()>c2->edgels.size()){
        return true;
    }
    return false;
 
}

static double get_continuity (const dbdet_edgel_chain *c1, const dbdet_edgel_chain *c2);
//: correct the CFG topology to produce a disjoint set


void dbdet_sel_base::merge_extreme_short_curve_frags()
{

for (unsigned i=0; i<edgemap_->edgels.size(); i++)
  {
    dbdet_edgel_chain *c1=0, *c2=0;
    dbdet_edgel* eA = edgemap_->edgels[i];

    int deg = curve_frag_graph_.pFrags[i].size()+ curve_frag_graph_.cFrags[i].size();	
    if (deg<2)
      continue; //nodes

    if (deg==2){ //degree 2 segments will trigger a splice
   
      //standard operation: extract them from the graph, reorder them, either merge or put them back

      //segments need to meet continuity criteria (simple one for now)
      if (curve_frag_graph_.pFrags[i].size()>1){
        dbdet_edgel_chain_list_iter fit = curve_frag_graph_.pFrags[i].begin();
        c1 =  (*fit); fit++;
        c2 =  (*fit);

        curve_frag_graph_.extract_fragment(c1);
        curve_frag_graph_.extract_fragment(c2);

        //reverse the sequence of edgels
        vcl_reverse(c2->edgels.begin(), c2->edgels.end());
        curve_frag_graph_.insert_fragment(c1);
        curve_frag_graph_.insert_fragment(c2);
      }
      else if (curve_frag_graph_.pFrags[i].size()==1){
        c1 =  curve_frag_graph_.pFrags[i].front();
        c2 =  curve_frag_graph_.cFrags[i].front();
    //for the closed contour case
      if(c1==c2)
        continue;
      }
      else {
        dbdet_edgel_chain_list_iter fit = curve_frag_graph_.cFrags[i].begin();
        c1 =  (*fit); fit++;
        c2 =  (*fit);

        //add the second one to the first one and delete it from the graph
        curve_frag_graph_.extract_fragment(c1);
        curve_frag_graph_.extract_fragment(c2);

        //reverse the sequence of edgels
        vcl_reverse(c1->edgels.begin(), c1->edgels.end());
        curve_frag_graph_.insert_fragment(c1);
        curve_frag_graph_.insert_fragment(c2);
      }

      if ((c1->edgels.size()<=5 || c2->edgels.size()<=5) && is_continue(c1,c2)){ //if two contours are all very short < 5 edges and are continuous, merge them anyway
        //merge the two contours
    	curve_frag_graph_.extract_fragment(c1);
        curve_frag_graph_.extract_fragment(c2);
        c1->append(c2->edgels);
        curve_frag_graph_.insert_fragment(c1);
    //when it makes a closed contour, just count as the child frag rather than parent frag
        if(c1->edgels.front()==c1->edgels.back())
        curve_frag_graph_.pFrags[c1->edgels.front()->id].remove(c1);
        delete c2;
      }
	  else if (c1->edgels.size()<=2 || c2->edgels.size()<=2)
	  {
		//merge the two contours
    	curve_frag_graph_.extract_fragment(c1);
        curve_frag_graph_.extract_fragment(c2);
        c1->append(c2->edgels);
        curve_frag_graph_.insert_fragment(c1);
    	//when it makes a closed contour, just count as the child frag rather than parent frag
        if(c1->edgels.front()==c1->edgels.back())
        curve_frag_graph_.pFrags[c1->edgels.front()->id].remove(c1);
        delete c2;		
	  }
     }
  }

    // use the filter to prun out local problems again
    regular_contour_filter();
}

void dbdet_sel_base::correct_CFG_topology()
{
  //D) Final T-junction type disambiguation can be done on the CFG
  // Basically, go over all the nodes of the CFG and operate on the ones with degree>2
  // also merge all segments that are adjacent

  // going over the edgemap instead so that a node is visited only once and so that I
  // don't have to deal with iterator issues

  for (unsigned i=0; i<edgemap_->edgels.size(); i++)
  {
    dbdet_edgel_chain *c1=0, *c2=0;
    dbdet_edgel* eA = edgemap_->edgels[i];

    int deg = curve_frag_graph_.pFrags[i].size()+ curve_frag_graph_.cFrags[i].size();	
    if (deg<2)
      continue; //nodes

    if (deg==2){ //degree 2 segments will trigger a splice
   
      //standard operation: extract them from the graph, reorder them, either merge or put them back

      //segments need to meet continuity criteria (simple one for now)
      if (curve_frag_graph_.pFrags[i].size()>1){
        dbdet_edgel_chain_list_iter fit = curve_frag_graph_.pFrags[i].begin();
        c1 =  (*fit); fit++;
        c2 =  (*fit);

        curve_frag_graph_.extract_fragment(c1);
        curve_frag_graph_.extract_fragment(c2);

        //reverse the sequence of edgels
        vcl_reverse(c2->edgels.begin(), c2->edgels.end());
        curve_frag_graph_.insert_fragment(c1);
        curve_frag_graph_.insert_fragment(c2);
      }
      else if (curve_frag_graph_.pFrags[i].size()==1){
        c1 =  curve_frag_graph_.pFrags[i].front();
        c2 =  curve_frag_graph_.cFrags[i].front();
    //for the closed contour case
    if(c1==c2)
        continue;
      }
      else {
        dbdet_edgel_chain_list_iter fit = curve_frag_graph_.cFrags[i].begin();
        c1 =  (*fit); fit++;
        c2 =  (*fit);

        //add the second one to the first one and delete it from the graph
        curve_frag_graph_.extract_fragment(c1);
        curve_frag_graph_.extract_fragment(c2);

        //reverse the sequence of edgels
        vcl_reverse(c1->edgels.begin(), c1->edgels.end());
        curve_frag_graph_.insert_fragment(c1);
        curve_frag_graph_.insert_fragment(c2);
      }

      if (is_continue(c1,c2)){ //if two contours are generally/local continue based on a frag
        //merge the two contours
    	curve_frag_graph_.extract_fragment(c1);
        curve_frag_graph_.extract_fragment(c2);
        c1->append(c2->edgels);
        curve_frag_graph_.insert_fragment(c1);
    //when it makes a closed contour, just count as the child frag rather than parent frag
        if(c1->edgels.front()==c1->edgels.back())
        curve_frag_graph_.pFrags[c1->edgels.front()->id].remove(c1);
        delete c2;
      }
     }
  }
/*
  // deal with junction with another loop after dealing with deg 2

  for (unsigned i=0; i<edgemap_->edgels.size(); i++)
  {
    dbdet_edgel_chain *c1=0, *c2=0;
    dbdet_edgel* eA = edgemap_->edgels[i];

    int deg = curve_frag_graph_.pFrags[i].size()+ curve_frag_graph_.cFrags[i].size();
    //degree 3 is a junction (T-junction or Y-junction)
    if (deg>2) //  Make length of contour as the first priority
    {
        //goal is to see if any two will produce smooth continuation
    dbdet_edgel_chain_list node_frags;
        if(curve_frag_graph_.pFrags[i].size()!=0){
        dbdet_edgel_chain_list_iter p_fit = curve_frag_graph_.pFrags[i].begin();
            for(;p_fit!=curve_frag_graph_.pFrags[i].end();p_fit++)
        {
            node_frags.push_back(*p_fit);
        }
    }

        if(curve_frag_graph_.cFrags[i].size()!=0){
        dbdet_edgel_chain_list_iter c_fit = curve_frag_graph_.cFrags[i].begin();
            for(;c_fit!=curve_frag_graph_.cFrags[i].end();c_fit++)
        {
            node_frags.push_back(*c_fit);
        }
    }
 
    node_frags.sort(is_longer); //sort all the pfrags and cfrags in length

        //compare each pair to decide merge or not
        dbdet_edgel_chain_list_iter fit_1=node_frags.begin();
    for (;fit_1!=--node_frags.end();){
        c1= *fit_1;
 
        if(c1->edgels.back()!= eA){ 
            curve_frag_graph_.extract_fragment(c1); 
            vcl_reverse(c1->edgels.begin(), c1->edgels.end());
            curve_frag_graph_.insert_fragment(c1);             
        }
        fit_1++;
        dbdet_edgel_chain_list_iter fit_2 = fit_1, max_fit = fit_1;
        double max_SM = 0;
        for (;fit_2!=node_frags.end();fit_2++){
            c2=*fit_2;
            if(c2->edgels.back()== eA){
                curve_frag_graph_.extract_fragment(c2);     
                vcl_reverse(c2->edgels.begin(), c2->edgels.end());
                curve_frag_graph_.insert_fragment(c2);             
            }

            double SM_0 = get_continuity(c1,c2);
            if(SM_0>max_SM){
                max_SM = SM_0;
                max_fit = fit_2;
            }
        }
        if(max_SM>=0.9){
            c2=*max_fit;
            curve_frag_graph_.extract_fragment(c1);
            curve_frag_graph_.extract_fragment(c2);
            c1->append(c2->edgels);
            curve_frag_graph_.insert_fragment(c1); 
            break;
        }
         
        //if(fit_2!=node_frags.end())
        //    break;
        }
    }
  }
*/
    // use the filter to prun out local problems again
    regular_contour_filter();
    vcl_cout<<"Finish correcting the CFG topology"<<vcl_endl;
}

//June, 2012:: Make decision whether to merge two contours or not (by Naman Kumar) 
static bool is_continue (const dbdet_edgel_chain *c1, const dbdet_edgel_chain *c2)
{
 	dbdet_edgel* e1;dbdet_edgel* e4;
	dbdet_edgel* e2 = c1->edgels.back(); // Last edge of c1
	dbdet_edgel* e3 = c2->edgels.front(); // Front edge of c2
	double dx1=0,dy1=0,dy2=0,dx2=0,s1=0,s2=0,s=0,SM_1=0;
	int j=0;
	for(int i=c1->edgels.size()-2;i>=0;i--)
	{
		e1 = c1->edgels[i];
    		j++;
    		if(j==c2->edgels.size()) break;
        	e4 = c2->edgels[j];
        	dx1 = e2->pt.x()-e1->pt.x();
        	dy1 = e2->pt.y()-e1->pt.y();
        	dx2 = e4->pt.x()-e3->pt.x();
        	dy2 = e4->pt.y()-e3->pt.y();
        	s1+=(e1)->strength; //strength1
        	s2+=(e4)->strength; //strength2
        	s=vcl_fabs((s1-s2)/j); // difference in strength
        	SM_1 = (dx1*dx2 + dy1*dy2)/vcl_sqrt(dx1*dx1+dy1*dy1)/vcl_sqrt(dx2*dx2+dy2*dy2);
//        	if(SM_1>=Theta_1 || s<=Strength_Diff1) return true;
//        	else if(SM_1<Theta_2 || s>Strength_Diff2) return false;

        	if(SM_1>=Theta_1) return true;
        	else if(SM_1<Theta_2) return false;

		else continue;
	}
    	return false;   
}

static double get_continuity (const dbdet_edgel_chain *c1, const dbdet_edgel_chain *c2)
{
        // using the median global continuity
        dbdet_edgel* e1=0;
        dbdet_edgel* e2=0;
        dbdet_edgel* e3=0;
        dbdet_edgel* e4=0;
        if(c1->edgels.size()>=5){ 
            e1 = c1->edgels[c1->edgels.size()-5];
            e2 = c1->edgels.back();
        }
        else{
            e1 = c1->edgels.front();
            e2 = c1->edgels.back();
        }
        if(c2->edgels.size()>=5){
            e3 = c2->edgels.front();
            e4 = c2->edgels[4];
        }
        else{
            e3 = c2->edgels.front();
            e4 = c2->edgels.back();
        }
        double dx1 = e2->pt.x()-e1->pt.x();
        double dy1 = e2->pt.y()-e1->pt.y();
        double dx2 = e4->pt.x()-e3->pt.x();
        double dy2 = e4->pt.y()-e3->pt.y();
        return (dx1*dx2 + dy1*dy2)/vcl_sqrt(dx1*dx1+dy1*dy1)/vcl_sqrt(dx2*dx2+dy2*dy2);
}

static bool share_same_ends(dbdet_edgel_chain *c1, dbdet_edgel_chain *c2)
{
    if((c1->edgels.front()==c2->edgels.front()&&c1->edgels.back()==c2->edgels.back())
        || (c1->edgels.front()==c2->edgels.back()&&c1->edgels.back()==c2->edgels.front()))
        return true;
    return false;
}



void dbdet_sel_base::regular_contour_filter(){
// first, deal with contours with length 3, which cause a lot of local problems
    dbdet_edgel_chain_list Size_3_chain_list;
    dbdet_edgel_chain_list_iter fit = curve_frag_graph_.frags.begin();
    for(;fit!=curve_frag_graph_.frags.end();fit++)
    {
        dbdet_edgel_chain *c1=*fit;
		// (a) for size 4 frags, git rid of the small closed triangle
		if(c1->edgels.size()==4 && (c1->edgels.front()==c1->edgels.back()))
		{
			c1->edgels.pop_back();
			curve_frag_graph_.pFrags[c1->edgels.back()->id].push_back(c1);
		}
		// push size 3 frags's pointer into a seperate list to same computation for next step
		if(c1->edgels.size()==3)
			Size_3_chain_list.push_back(c1);
    }


    dbdet_edgel_chain_list_iter fit_1 = Size_3_chain_list.begin();
    while(fit_1!=Size_3_chain_list.end())
    {
        // (b) change 3 edgels sharp path to be a direct line link
    dbdet_edgel_chain *c1=*fit_1;
    double dx1 = c1->edgels[1]->pt.x() - c1->edgels.front()->pt.x();
    double dy1 = c1->edgels[1]->pt.y() - c1->edgels.front()->pt.y();
    double dx2 = c1->edgels.back()->pt.x() - c1->edgels[1]->pt.x();
    double dy2 = c1->edgels.back()->pt.y() - c1->edgels[1]->pt.y();
    double SM = (dx1*dx2 + dy1*dy2)/vcl_sqrt(dx1*dx1+dy1*dy1)/vcl_sqrt(dx2*dx2+dy2*dy2);
    double length = vcl_sqrt(dx1*dx1+dy1*dy1) + vcl_sqrt(dx2*dx2+dy2*dy2);
    if(length>10)// only consider local problems
        continue;
    if(SM<=0)
    {
        dbdet_edgel_list_iter eit = c1->edgels.begin();
        eit++;
        c1->edgels.erase(eit);
    }
    // (c) if two frags share same end nodes, change them into one direct line link
    dbdet_edgel_chain_list_iter fit_2 = ++fit_1;
    while(fit_2!=Size_3_chain_list.end())
    {
        dbdet_edgel_chain *c2=*fit_2;
        if(share_same_ends(c1, c2))
        {
            // check if c1 is modifid in (b), if no, change it to a direct line link
            if(c1->edgels.size()==3)
            {
                dbdet_edgel_list_iter eit = c1->edgels.begin();
                eit++;
                c1->edgels.erase(eit);
            }
            // remove c2 from CFG
            curve_frag_graph_.extract_fragment(c2);
            // only deal with one frag has same ends, if there are more, leave for iterations afterward
            break;
        }
        fit_2++;
    }
    }
}


//New construction of Hypothesis Tree by Naman Kumar
void dbdet_sel_base::Construct_Hypothesis_Tree()
{



	regular_contour_filter(); // Filtering step
	int n1=0;
	double d=0,dis=0,distance=0;
	vcl_vector<dbdet_edgel*> new_chain0,new_chain2,new_chain3,new_chain6,new_chain33;
	dbdet_edgel_chain* new_chain1=new dbdet_edgel_chain();dbdet_edgel_chain* new_chain4=new dbdet_edgel_chain();
	dbdet_edgel_chain* test1=new dbdet_edgel_chain();dbdet_edgel_chain *chains=new dbdet_edgel_chain();
	dbdet_edgel_chain* new_chain44=new dbdet_edgel_chain();
	double gap_thres=gap_;
	vcl_cout << "Construction of Hypothesis Tree is in Progress!! " << vcl_endl;
	//Calculating number of edges which are having degree 1   	
	for (int i=0; i<edgemap_->edgels.size(); i++)
	{
        	dbdet_edgel* eA1 = edgemap_->edgels[i];
        	new_chain0.push_back(eA1);
        	if ((curve_frag_graph_.pFrags[eA1->id].size() + curve_frag_graph_.cFrags[eA1->id].size()) ==1) new_chain1->edgels.push_back(eA1);
        	else new_chain2.push_back(eA1);
	}

	//Calculating number of edges which are part of the contours and which are unused
	dbdet_edgel_chain_list_iter fit = curve_frag_graph_.frags.begin();
	for(;fit!=curve_frag_graph_.frags.end();fit++)
	{
        	dbdet_edgel_chain *test1=*fit;
        	for(int d=0;d<test1->edgels.size();d++) {new_chain3.push_back(test1->edgels[d]);new_chain33.push_back(test1->edgels[d]);} 
	}
	
	for(int i=0;i<new_chain0.size();i++)
	{
        	for(int j=0;j<new_chain3.size();j++)
        	{
        		if(new_chain0[i]!=new_chain3[j]) continue;
                	else {n1=5; break;}
        	}
     		if(n1==0) {new_chain4->edgels.push_back(new_chain0[i]);}
     		else {n1=0; continue;}
	}
	//Constucting the tree from end of an unambiguous chain and extending it till the end of edge chain
	double cost1=gap_,cost2=10.0,cost3=gap_,d1=0.0,d2=0.0,d3=0,dx=0.0,dy=0.0,cost=1000.0,costc=0.0;int m1=0,m2=0,m3=0,m4=0,m5=0,m7=0,m8=0,m9=0;
	dbdet_edgel* ce=0;dbdet_edgel* pe=0;dbdet_edgel* ed=0;dbdet_edgel* imp=0;dbdet_edgel* im=0;
	dbdet_edgel_chain *c11=new dbdet_edgel_chain();dbdet_edgel_chain* xx=new dbdet_edgel_chain();dbdet_edgel_chain* end=new dbdet_edgel_chain();
	while(new_chain1->edgels.size()>0)
    	{
        	a: ;
        	if(new_chain1->edgels.size()==0) break; 
        	ed=new_chain1->edgels[0];     
        	new_chain1->edgels.pop_front();
        	for(int z=0;z<end->edgels.size();z++){if(ed==end->edgels[z]) goto a; else continue;}
        	dbdet_edgel_chain* new_chain5 = new dbdet_edgel_chain();
        	new_chain5->edgels.push_back(ed);xx->edgels.push_back(ed);
        	m4=0;m7=0;
		double dis=0, distance=0;
        	//Forming the tree from the edge
        	while(1)
        	{
                	dbdet_edgel_list_iter eit1=new_chain4->edgels.begin(); dbdet_edgel_list_iter eit2=new_chain4->edgels.begin();
                 	if(m4==0)
                 	{
                    		if(curve_frag_graph_.cFrags[ed->id].size()==1) // if number of child fragments is 1
                    		{	dbdet_edgel_chain_list_iter ccit = curve_frag_graph_.cFrags[ed->id].begin();
                        	 	ce = (*ccit)->edgels[1];c11=*ccit;pe=ce;m7=1;			
					for (int j=1;j<c11->edgels.size();j++) 
					{dis=vgl_distance(c11->edgels[j]->pt,c11->edgels[j-1]->pt);if(dis>distance) distance=dis;}
			     		distance=distance + 0.25;if(distance <=1) gap_=1; else if(distance <gap_thres) gap_=distance; else gap_=gap_thres; 
                    		}
                    		else if(curve_frag_graph_.pFrags[ed->id].size()==1) // if number of parent fragments is 1
                   		{	dbdet_edgel_chain_list_iter pcit = curve_frag_graph_.pFrags[ed->id].begin();
                        		ce = (*pcit)->edgels[(*pcit)->edgels.size()-2];c11=*pcit;pe=ce;m7=2;		
					for (int j=1;j<c11->edgels.size();j++) 
					{dis=vgl_distance(c11->edgels[j]->pt,c11->edgels[j-1]->pt);if(dis>distance) distance=dis;}
			     		distance=distance + 0.25;if(distance <=1) gap_=1; else if(distance <gap_thres) gap_=distance; else gap_=gap_thres;	
                    		}m4=1;
                 	} 
                	// Used later
                	if(m7==2){m8=c11->edgels.size()-5;if(m8<0) m8=0;m9=c11->edgels.size();}
                	else if(m7==1) {m8=0;m9=5; if(m9>c11->edgels.size())m9=c11->edgels.size();}
		
			// Finding the closest unused edge
                	costc=0.0;cost=10000.0;cost1=gap_;
                	for(int j=0;j<new_chain4->edgels.size(); j++)
                	{ 
                        	d1= vgl_distance(ed->pt,new_chain4->edgels[j]->pt);
                        	//Checking Localization, Orientation,etc..
                        	if(d1<cost1)
                        	{
                            		vcl_vector<dbdet_edgel*> dummy_chain;
					dbdet_edgel_chain* edgel_chain = new dbdet_edgel_chain();
                                  	for(int i=0;i<new_chain5->edgels.size();i++)edgel_chain->edgels.push_back(new_chain5->edgels[i]);
                                  	edgel_chain->edgels.push_back(new_chain4->edgels[j]);
                                  	vcl_vector<dbdet_edgel*> chain(edgel_chain->edgels.begin(),edgel_chain->edgels.end());
                                  	costc = compute_path_metric2(dummy_chain, chain, dummy_chain);
                                  	if(costc<cost)
                                   	{
                         			double d8=vgl_distance(new_chain4->edgels[j]->pt,ce->pt);
                                      		double d9=vgl_distance(new_chain4->edgels[j]->pt,ed->pt);
                                        	double d0=vgl_distance(new_chain4->edgels[j]->pt,pe->pt);
                                        	double dx1 = ce->pt.x() - ed->pt.x();
                                        	double dy1 = ce->pt.y() - ed->pt.y();
                                        	double dx2 = ed->pt.x() - new_chain4->edgels[j]->pt.x();
                                        	double dy2 = ed->pt.y() - new_chain4->edgels[j]->pt.y();
						double angle=((dx1*dx2 + dy1*dy2)/vcl_sqrt(dx1*dx1+dy1*dy1)/vcl_sqrt(dx2*dx2+dy2*dy2));
						if(d0<d9) {++eit1;continue;}
                                        	if(d8<d9 || angle<0){++eit1;continue;}
                                        	imp=new_chain4->edgels[j];
                                        	cost=costc;
                                        	m1=1;
                                        	eit2=eit1;
                                   	}  				          
                         	}
                         	else if(d1<cost2 && d1>1){cost2=d1;}
                         	++eit1;     
                  	} 
          
                  	m3=0;m5=0;cost3=gap_;
                  	// Finding the closest edge which is part of a fragment
                  	for(int t=0;t<new_chain3.size(); t++)
                  	{
                  		if(new_chain3[t]==ed || new_chain3[t]==ce) continue;
                    		d1= vgl_distance(ed->pt,new_chain3[t]->pt);
                        	if(d1<=cost3)
                        	{
					//Dont consider the previous 5 edges present in parent/child fragment of the starting edge
                                	for(int c=m8;c<m9;c++) {if(new_chain3[t]==c11->edgels[c]) goto z; else continue;}
                                	//Dont use the edge which is part of the same tree again 
                                	for(int c=0;c<new_chain5->edgels.size();c++) {if(new_chain3[t]==new_chain5->edgels[c]) goto z; else continue;} 
                                	im=new_chain3[t];     
                                	cost3=d1;m5=1;
                                	dx=vgl_distance(im->pt,ce->pt);
                                	dy=vgl_distance(im->pt,ed->pt);
                        	}     
                   		z: ;
                  	}
              		if(dx>dy && m5==1){m3=5;m1=1;imp=im;}
                  	if(m1==1)
                  	{
                   		m2=1;m1=0;cost1=gap_;
                    		ce=ed;ed=imp; xx->edgels.push_back(imp);new_chain5->edgels.push_back(imp);
                    		if(m3==0){new_chain4->edgels.erase(eit2);new_chain3.push_back(imp);new_chain33.push_back(imp);}
                    		if(m3!=0){m3=0;break;}
                  	}
                  	else if(cost2>1) break;
        	}
        	//No double contours within the same 2 end points.
        	if(m2==1) {if(c11->edgels.front()==c11->edgels.back()) continue;} 
        	//Add the tree   
        	new_chain5->temp = true;
        	curve_frag_graph_.CFTG.insert_fragment(new_chain5);
        	end->edgels.push_back(new_chain5->edgels.back());
   
	}
	dbdet_edgel* edge1=0;dbdet_edgel* edge2=0;dbdet_edgel_chain *chain1=new dbdet_edgel_chain();vcl_list<dbdet_CFTG_link*> GD_list;
	double p1=1.0;int p2=0,p3=0,p16=0;
	dbdet_CFTG_link_list_iter l_it = curve_frag_graph_.CFTG.Links.begin();      
    for (; l_it != curve_frag_graph_.CFTG.Links.end(); l_it++) GD_list.push_back(*l_it);
    //while size is greater than 0
  	while (GD_list.size()>0)
  	{
    double distance1=0,distance2=0;
		dbdet_CFTG_link* cur_Link = GD_list.front();GD_list.pop_front();
    dbdet_edgel_chain_list_iter f_it = cur_Link->cCFs.begin();dbdet_edgel_chain* new_chain5=(*f_it);dbdet_edgel* edge3=new_chain5->edgels.front();
		gap_=1;
		// iterating through all the edgels of chain5
    for(int i=0;i<new_chain5->edgels.size();i++)
		{
			dbdet_edgel_chain *new_chain6a=new dbdet_edgel_chain();
			p1=gap_;
			dbdet_edgel_list_iter eit5=new_chain4->edgels.begin(); dbdet_edgel_list_iter eit6=new_chain4->edgels.begin();
			for(int j=0;j<new_chain4->edgels.size();j++)
			{
				double p4=vgl_distance(new_chain5->edgels[i]->pt,new_chain4->edgels[j]->pt);				
				if(p4<p1){edge1=new_chain4->edgels[j];p1=p4;p2=1;eit6=eit5;++eit5;} // if distance between edgels is less than gap
				else {++eit5;continue;}	
			}
			if(p2==1)
			{
				new_chain4->edgels.erase(eit6); //remove edgel
				dbdet_edgel* edge4=0;p2=0;double p5=gap_;
				for(int b=0;b<new_chain5->edgels.size();b++)
				{
					double p6=vgl_distance(edge1->pt,new_chain5->edgels[b]->pt);
					if(p6<p5){edge2=new_chain5->edgels[b];p5=p6;} 
					else continue;
				}
				new_chain6a=new dbdet_edgel_chain();new_chain6a->edgels.push_back(edge2);new_chain6a->edgels.push_back(edge1);
				int p7=0,p8=0,p9=0;
                		if(curve_frag_graph_.cFrags[edge3->id].size()>=1) // if size of child fragment is greater than 1
                                {dbdet_edgel_chain_list_iter ccit = curve_frag_graph_.cFrags[edge3->id].begin();chain1=(*ccit);p7=1;}
                            	else if(curve_frag_graph_.pFrags[edge3->id].size()>=1) // if size of parent fragment is greater than 1
                             	{dbdet_edgel_chain_list_iter pcit = curve_frag_graph_.pFrags[edge3->id].begin();chain1=(*pcit);p7=2;}
                		if(p7==2){p8=chain1->edgels.size()-5;if(p8<0) p8=0;p9=chain1->edgels.size();}
                        	else if(p7==1) {p8=0;p9=5; if(p9>chain1->edgels.size())p9=chain1->edgels.size();}
				while(1)
				{
					double p10=gap_,p11=0.0;
					dbdet_edgel_list_iter eit3=new_chain4->edgels.begin();dbdet_edgel_list_iter eit4=new_chain4->edgels.begin();	
          // iterating through all the edgels of the chain
          for(int a=0;a<new_chain4->edgels.size();a++)
					{
						p11= vgl_distance(edge1->pt,new_chain4->edgels[a]->pt);
						if(p11<p10)
						{
							double d8=vgl_distance(new_chain4->edgels[a]->pt,edge2->pt);
              double d9=vgl_distance(new_chain4->edgels[a]->pt,edge1->pt);
              double dx1 = edge1->pt.x() - edge2->pt.x(), dy1 = edge1->pt.y() - edge2->pt.y();
              double dx2=new_chain4->edgels[a]->pt.x()-edge1->pt.x();
							double dy2=new_chain4->edgels[a]->pt.y()-edge1->pt.y();
              if(d8<d9 || ((dx1*dx2 + dy1*dy2)/vcl_sqrt(dx1*dx1+dy1*dy1)/vcl_sqrt(dx2*dx2+dy2*dy2))<0.4) 
							{++eit3;continue;}
							p10=p11;edge4=new_chain4->edgels[a];eit4=eit3;p2=1;
						}
						++eit3;
					}
					double p12=gap_,p13=0,p14=0.0,p15=0.0;p16=0;dbdet_edgel* edge5=0;
					for(int t=0;t<new_chain3.size(); t++)
                  			{
						if(new_chain3[t]==edge1) continue;                    			       		
						p13= vgl_distance(edge1->pt,new_chain3[t]->pt);
            if(p13<=p12)
              {
              for(int c=p8;c<p9;c++) {if(new_chain3[t]==chain1->edgels[c])goto jump;else continue;} 						      		
              for(int a=0;a<new_chain5->edgels.size();a++)
							{if(new_chain3[t]==new_chain5->edgels[a])goto jump; else continue;}
					     for(int b=0;b<new_chain6a->edgels.size();b++)
							{if(new_chain3[t]==new_chain6a->edgels[b])goto jump;else continue;}
              edge5=new_chain3[t];p12=p13;p14=vgl_distance(edge5->pt,edge2->pt);
							p15=vgl_distance(edge5->pt,edge1->pt);
					    p16=1;
              }     
                   			       	jump: ;
                  			}
					if(p14>p15 && p16==1) {edge4=edge5;p2=1;}					
					if(p2==1)
					{
						new_chain6a->edgels.push_back(edge4);p3=1;p2=0;
						if(p16==0){new_chain4->edgels.erase(eit4);new_chain3.push_back(edge4);edge2=edge1;edge1=edge4;}else break;
					}
					else break;		
				}
			        double p17=0,p18=0,p21=0;
				if(p3==1 && new_chain6a->edgels.size()>5) // minimum size should be 5
				{
					dbdet_edgel* edge6=new_chain6a->edgels[new_chain6a->edgels.size()/2];				
					for(int i=0;i<new_chain33.size();i++)
					{
						p18=vgl_distance(edge6->pt,new_chain33[i]->pt);if(p18<1) p21=10; // if distance is less than 1
						if(p21==10){p17=1;break;}
					}	
					if(p17==0){new_chain6a->temp = true;curve_frag_graph_.CFTG.insert_fragment(new_chain6a);} // insert the fragment
					p3=0;p2=0;p16=0,p3=0;
				}
			}
		}
	}

// Add by Yuliang, a indicator shows edges participating in unambiguous frags and hypothesis trees
	vcl_cout << "counting in participate edges" << vcl_endl;

	// count in hypothesis tree edges
	dbdet_CFTG_link_list_iter it_8 = curve_frag_graph_.CFTG.Links.begin();

	for (; it_8 != curve_frag_graph_.CFTG.Links.end(); it_8++)
	{
		dbdet_edgel_chain_list cur_cCFs = (*it_8)->cCFs;
		
		dbdet_edgel_chain_list_iter it_88 = cur_cCFs.begin();
		for(; it_88 != cur_cCFs.end(); it_88++)
		{
			dbdet_edgel_list_iter it_888 = (*it_88)->edgels.begin();
			for(; it_888 != (*it_88)->edgels.end(); it_888++)
				curve_frag_graph_.participate_edge_id.insert((*it_888)->id);
		}
	}



/*	std::set<int>::iterator it = curve_frag_graph_.participate_edge_id.begin();
	for (; it!=curve_frag_graph_.participate_edge_id.end(); it++)
		vcl_cout << (*it) << " ";
	vcl_cout << vcl_endl;
*/

	vcl_cout <<"participate edge count: " <<curve_frag_graph_.participate_edge_id.size() << vcl_endl;
	
	vcl_cout << "Hypothesis Tree Constructed!!" << vcl_endl;	
}

// New Disambiguation Process by Naman Kumar
void dbdet_sel_base::Disambiguation()
{
	vcl_cout << "Disambiguating the Hypothesis Tree!!" << vcl_endl;
	dbdet_CFTG_link_list_iter l_it = curve_frag_graph_.CFTG.Links.begin();
	for (; l_it != curve_frag_graph_.CFTG.Links.end(); l_it++)
	{
        	double cost=0.0;
        	int deg_S = curve_frag_graph_.CFTG.cLinks[(*l_it)->eS->id].size() + curve_frag_graph_.CFTG.pLinks[(*l_it)->eS->id].size();
        	dbdet_CFTG_link* cur_Link = (*l_it);
		//Calculating the Cost        	
		vcl_vector<dbdet_edgel*> dummy_chain;
        	dbdet_edgel_chain* edgel_chain = cur_Link->cCFs.front();
        	vcl_vector<dbdet_edgel*> chain(edgel_chain->edgels.begin(),edgel_chain->edgels.end());
	 	cost = compute_path_metric2(dummy_chain, chain, dummy_chain);
        	//Degree = 1
        	if(deg_S==1) {curve_frag_graph_.insert_fragment((*l_it)->cCFs.front()); continue;}
        	//Degree > 1
        	// To fill in the small gaps in closed contours
		if((curve_frag_graph_.pFrags[edgel_chain->edgels.front()->id].size()+curve_frag_graph_.cFrags[edgel_chain->edgels.front()->id].size())==1    			&&(curve_frag_graph_.pFrags[edgel_chain->edgels.back()->id].size()+curve_frag_graph_.cFrags[edgel_chain->edgels.back()->id].size())==1 			&& edgel_chain->edgels.size()==2) curve_frag_graph_.insert_fragment(edgel_chain);
        	if(cost<1.0 && edgel_chain->edgels.size()>2) curve_frag_graph_.insert_fragment(edgel_chain); 
	}
	//clear the graph
	curve_frag_graph_.CFTG.clear();
	curve_frag_graph_.CFTG.resize(edgemap_->edgels.size());
}

// By Naman Kumar: a minor function just to prune some extra small part of contours
void dbdet_sel_base::Post_Process()
{

	dbdet_edgel_chain* new_chain= new dbdet_edgel_chain();
  // iterating through all the contours and deleting those that have length less than a certain threshold
	for (int i=0; i<edgemap_->edgels.size(); i++)
    	{
        	dbdet_edgel* eA1 = edgemap_->edgels[i];
            	if ((curve_frag_graph_.pFrags[eA1->id].size() + curve_frag_graph_.cFrags[eA1->id].size()) ==1) new_chain->edgels.push_back(eA1);
    	}   
  for(int j=0;j<new_chain->edgels.size();j++)
    	{
        	dbdet_edgel* edge=new_chain->edgels[j]; dbdet_edgel* edge2=0;dbdet_edgel_chain* chain= new dbdet_edgel_chain();dbdet_edgel* edge3=0;
        	int n=0,number=0,num=0,diff=0;dbdet_edgel* edge4=0;
        	if(curve_frag_graph_.cFrags[edge->id].size()==1) // if size of child fragment is 1
            	{
            		n=1;dbdet_edgel_chain_list_iter ccit = curve_frag_graph_.cFrags[edge->id].begin();chain=*ccit;
            		edge2 = chain->edgels[1];if(chain->edgels.size()>2){number=1;edge3=chain->edgels[2];}diff=1;
		          }
            	else if(curve_frag_graph_.pFrags[edge->id].size()==1) // if size of parent fragment is 1
            	{
            		n=1;dbdet_edgel_chain_list_iter pcit = curve_frag_graph_.pFrags[edge->id].begin();chain=*pcit;
            		edge2 = chain->edgels[chain->edgels.size()-2];
            		if(chain->edgels.size()>2){number=1;edge3=chain->edgels[chain->edgels.size()-3];}diff=2;
		          } 

        	dbdet_edgel_chain_list_iter fit = curve_frag_graph_.frags.begin();
          // going through all the curve fragments of the graph
        	for(;fit!=curve_frag_graph_.frags.end();fit++)
            	{dbdet_edgel_chain *test1=*fit;if(test1==chain)continue;for(int d=0;d<test1->edgels.size();d++){if(edge==test1->edgels[d]) goto end;}}
        	if(n==1 && chain->edgels.size()>1 && ((curve_frag_graph_.pFrags[edge2->id].size() + curve_frag_graph_.cFrags[edge2->id].size()) >=1))
        	{
            		curve_frag_graph_.extract_fragment(chain); // extract the fragment
            		if(diff==1) chain->edgels.pop_front();else if(diff==2) chain->edgels.pop_back();curve_frag_graph_.insert_fragment(chain);

        	}
        	else if(number==1 && chain->edgels.size()>1 && ((curve_frag_graph_.pFrags[edge3->id].size()+curve_frag_graph_.cFrags[edge3->id].size()) >=1))
        	{           
            		curve_frag_graph_.extract_fragment(chain); // extract the fragment
                	if(diff==1) {chain->edgels.pop_front();chain->edgels.pop_front();}
            		else if(diff==2){chain->edgels.pop_back();chain->edgels.pop_back();}curve_frag_graph_.insert_fragment(chain);
        	}
       		end: ;
    	}
} 
