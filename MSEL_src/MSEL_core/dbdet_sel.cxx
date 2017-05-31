#include "dbdet_sel.h"

#include <vcl_iostream.h>
#include <vcl_fstream.h>
#include <vcl_cassert.h>
#include <vcl_deque.h>
#include <vcl_algorithm.h>

//: form an edgel pair (ref_e->e2) or (e2->ref_e)
template <class curve_model>
void dbdet_sel<curve_model>::form_an_edgel_pair(dbdet_edgel* ref_e, dbdet_edgel* e2)
{
  //form a hypothetical curve model from this pair
  bool ref_first;
  curve_model* cm = form_a_hypothesis(ref_e, e2, ref_first);
  
  //if model is valid, form an edgel pair curvelet
  if (cm){
    dbdet_curvelet* c1 = new dbdet_curvelet(ref_e, cm);
    if (ref_first){
      c1->push_back(ref_e);
      c1->push_back(e2);
    }
    else {
      c1->push_back(e2);
      c1->push_back(ref_e);
    }
    
    //make links to the edgel-pair from the ref edgel
    curvelet_map_.add_curvelet(c1);
  } 
}

////: form an edgel triplet from two pairs
//void dbdet_sel::form_an_edgel_triplet(dbdet_curvelet* p1, dbdet_curvelet* p2)
//{
//  // combine them using syntactic rules
//
//  // RULE 1. AX + AY = AXY or AYX
//  if (p1->edgel_chain[0] == p2->edgel_chain[0])
//  { 
//    dbdet_edgel* eA = p1->edgel_chain[0];//common edgel
//    dbdet_edgel* eX = p1->edgel_chain[1];
//    dbdet_edgel* eY = p2->edgel_chain[1];
//
//    //compute intersection of curve bundles for each perturbation of this edgel
//    dbdet_curve_bundle* eAeb = new dbdet_curve_bundle(eA, p1->curve_bundles[0], p2->curve_bundles[0]);
//
//    //if any of the perturbations produced a valid bundle, form an edgel triplet
//    if (eAeb->any_valid_bundles()){
//      dbdet_curvelet* trip = new dbdet_curvelet();
//
//      // decide between AXY or AYX using distance
//      if (vgl_distance(eA->pt, eX->pt) < vgl_distance(eA->pt, eY->pt)){
//        trip->push_back(eA, eAeb);
//        trip->push_back(eX, 0);      //no curve bundles for this edgel
//        trip->push_back(eY, 0);      //no curve bundles for this edgel
//      } 
//      else {
//        trip->push_back(eA, eAeb);
//        trip->push_back(eY, 0);      //no curve bundles for this edgel
//        trip->push_back(eX, 0);      //no curve bundles for this edgel
//      }
//
//      //record this triplet
//      add_a_triplet(trip, 0, p1, p2);
//    }
//  }
//
//  // RULE 2. AX + YA = YAX (this is one of the primary rules)
//  else if (p1->edgel_chain[0] == p2->edgel_chain[1])
//  { 
//    dbdet_edgel* eX = p1->edgel_chain[1];
//    dbdet_edgel* eA = p1->edgel_chain[0];//common edgel
//    dbdet_edgel* eY = p2->edgel_chain[0];
//
//    //compute intersection of curve bundles for each perturbation of this edgel
//    dbdet_curve_bundle* eAeb = new dbdet_curve_bundle(eA, p1->curve_bundles[0], p2->curve_bundles[1]);
//
//    //if any of the perturbations produced a valid bundle, form an edgel triplet
//    if (eAeb->any_valid_bundles()){
//      dbdet_curvelet* trip = new dbdet_curvelet();
//
//      trip->push_back(eY, 0);      //no curve bundles for this edgel
//      trip->push_back(eA, eAeb);
//      trip->push_back(eX, 0);      //no curve bundles for this edgel
//
//      //record this triplet
//      add_a_triplet(trip, 1, p2, p1);
//    }
//  }
//
//  // RULE 3. XA + AY = XAY (this is one of the primary rules)
//  else if (p1->edgel_chain[1] == p2->edgel_chain[0])
//  { 
//    dbdet_edgel* eX = p1->edgel_chain[0];
//    dbdet_edgel* eA = p1->edgel_chain[1];//common edgel
//    dbdet_edgel* eY = p2->edgel_chain[1];
//
//    //compute intersection of curve bundles for each perturbation of this edgel
//    dbdet_curve_bundle* eAeb = new dbdet_curve_bundle(eA, p1->curve_bundles[1], p2->curve_bundles[0]);
//
//    //if any of the perturbations produced a valid bundle, form an edgel triplet
//    if (eAeb->any_valid_bundles()){
//      dbdet_curvelet* trip = new dbdet_curvelet();
//
//      trip->push_back(eX, 0);      //no curve bundles for this edgel
//      trip->push_back(eA, eAeb);
//      trip->push_back(eY, 0);      //no curve bundles for this edgel
//
//      //record this triplet
//      add_a_triplet(trip, 1, p1, p2);
//    }
//  }
//
//  // RULE 4. XA + YA = XYA or YXA
//  else if (p1->edgel_chain[1] == p2->edgel_chain[1])
//  { 
//    dbdet_edgel* eA = p1->edgel_chain[1];//common edgel
//    dbdet_edgel* eX = p1->edgel_chain[0];
//    dbdet_edgel* eY = p2->edgel_chain[0];
//
//    //compute intersection of curve bundles for each perturbation of this edgel
//    dbdet_curve_bundle* eAeb = new dbdet_curve_bundle(eA, p1->curve_bundles[1], p2->curve_bundles[1]);
//
//    //if any of the perturbations produced a valid bundle, form an edgel triplet
//    if (eAeb->any_valid_bundles()){
//      dbdet_curvelet* trip = new dbdet_curvelet();
//
//      // decide between XYA or YXA using distance
//      if (vgl_distance(eA->pt, eX->pt) < vgl_distance(eA->pt, eY->pt)){
//        trip->push_back(eY, 0);      //no curve bundles for this edgel
//        trip->push_back(eX, 0);      //no curve bundles for this edgel
//        trip->push_back(eA, eAeb);
//      } 
//      else {
//        trip->push_back(eX, 0);      //no curve bundles for this edgel
//        trip->push_back(eY, 0);      //no curve bundles for this edgel
//        trip->push_back(eA, eAeb);
//      }
//
//      //record this triplet
//      add_a_triplet(trip, 2, p1, p2);
//    }
//  }
//}
//
////: add a newly formed edgel triplet to all relevant lists
//void dbdet_sel::add_a_triplet(dbdet_curvelet* trip, unsigned pos, 
//                               dbdet_curvelet* p1, dbdet_curvelet* p2)
//{
//  //need to do a triplet existence check first because this triplet might already 
//  //have formed due to some other rule
//  dbdet_curvelet* p1_prev_trip = p1->exists(trip);
//  dbdet_curvelet* p2_prev_trip = p2->exists(trip);
//
//  //if it exists in p1
//  if (p1_prev_trip)
//  {
//    //if it already exists check to see if the common edgel already has a curve bundle
//    if (!p1_prev_trip->curve_bundles[pos]){
//      p1_prev_trip->curve_bundles[pos] = trip->curve_bundles[pos];
//      trip->curve_bundles[pos] = 0;
//    }
//    delete trip; //this triplet is not required anymore
//
//    //if p2 is not linked to this triplet, make a link
//    if (!p2_prev_trip)
//      p2->add_larger_curvelet(p1_prev_trip);
//  }
//  else if (p2_prev_trip){ //if it exists in p2
//    //if it already exists check to see if the common edgel already has a curve bundle
//    if (!p2_prev_trip->curve_bundles[pos]){
//      p2_prev_trip->curve_bundles[pos] = trip->curve_bundles[pos];
//      trip->curve_bundles[pos] = 0;
//    }
//    delete trip; //this is not required anymore
//
//    //if p1 is not linked to this triplet, make a link
//    if (!p1_prev_trip)
//      p1->add_larger_curvelet(p2_prev_trip);
//  }
//  else { 
//    //brand new triplet
//
//    //add links from the edgels to this triplet
//    trip->add_links_from_participating_edgels();
//
//    //also make links from the participating pairs
//    p1->add_larger_curvelet(trip);
//    p2->add_larger_curvelet(trip);
//
//    //and finally, add it to the edge triplets list
//    trips.push_back(trip);
//  }
//}
//
////: remove redundant pairs after triplets have formed
//void dbdet_sel::remove_redundant_pairs()
//{
//  // if ABC exists, then AC is redundant. 
//  // so remove AC and all groups that AC formed.
//
//  vcl_set<vcl_pair<dbdet_edgel*, dbdet_edgel*> > redundant_pairs;
//
//  //go over the list of edgel triplets and find redundant pairs
//  vcl_list<dbdet_curvelet* >::iterator t_it = trips.begin();
//  for (; t_it != trips.end(); t_it++)
//  {
//    dbdet_curvelet* trip = (*t_it);
//
//    //record the pair AC
//    redundant_pairs.insert(
//      vcl_pair<dbdet_edgel*, dbdet_edgel*>(trip->edgel_chain[0], trip->edgel_chain[2]));
//
//  }
//
//  //now go over the redundant pairs list and check if any of them exist
//  vcl_set<vcl_pair<dbdet_edgel*, dbdet_edgel*> >::iterator rp_it = redundant_pairs.begin();
//  for (; rp_it != redundant_pairs.end(); rp_it++)
//  {
//    dbdet_edgel* eA = rp_it->first;
//    dbdet_edgel* eC = rp_it->second;
//
//    //find the current pair in the edgel's pairs list
//    dbdet_curvelet* pair = eA->find_pair(eA, eC);
//    if (pair){
//      //if this exists, go over all the triplets that it formed
//      curvelet_list_iter c_it = pair->larger_curvelets.begin();
//      for (; c_it != pair->larger_curvelets.end(); c_it++){
//        dbdet_curvelet* t1 = (*c_it);
//
//        //if any of these are of the form XAC or ACY, delete them
//        if ((t1->edgel_chain[0]==eA && t1->edgel_chain[1]==eC) ||
//            (t1->edgel_chain[1]==eA && t1->edgel_chain[2]==eC))
//        {
//          //remove links from the edgels
//          t1->edgel_chain[0]->remove_curvelet(t1);
//          t1->edgel_chain[1]->remove_curvelet(t1);
//          t1->edgel_chain[2]->remove_curvelet(t1);
//          //also remove links from other pairs
//
//          delete t1;
//        }
//      }
//      //finally delete the pair (AC)
//      eA->remove_curvelet(pair); //remove links from the edgels first
//      eC->remove_curvelet(pair);
//      delete pair;
//    }
//  }
//
//}
//
////: form an edgel quad from two triplets
//void dbdet_sel::form_an_edgel_quad(dbdet_curvelet* t1, dbdet_curvelet* t2)
//{
//  // combine them using syntactic rules
//
//  // RULE 1. ABX + ABY = ABXY or ABYX
//  if (t1->edgel_chain[0] == t2->edgel_chain[0] && 
//      t1->edgel_chain[1] == t2->edgel_chain[1])
//  { 
//    dbdet_edgel* eA = t1->edgel_chain[0];//common edgel
//    dbdet_edgel* eB = t1->edgel_chain[1];//common edgel
//    dbdet_edgel* eX = t1->edgel_chain[2];
//    dbdet_edgel* eY = t2->edgel_chain[2];
//
//    //compute intersection of curve bundles for each perturbation of the common edgels
//    dbdet_curve_bundle *eAcb=0, *eBcb=0;
//    if (t1->curve_bundles[0] && t2->curve_bundles[0])
//      eAcb = new dbdet_curve_bundle(eA, t1->curve_bundles[0], 
//                                        t2->curve_bundles[0] );
//
//    if (t1->curve_bundles[1] && t2->curve_bundles[1])
//      eBcb = new dbdet_curve_bundle(eB, t1->curve_bundles[1], 
//                                        t2->curve_bundles[1] );
//
//    //The participating curve bundles at at least one of the common edgels have to be defined in order to form a quad
//    if (!eAcb && !eBcb)
//      return;
//
//    //if any of the perturbations produced a valid bundle, form an edgel quad
//    if ((eAcb && eAcb->any_valid_bundles()) || 
//        (eBcb && eBcb->any_valid_bundles()))
//    {
//      dbdet_curvelet* quad = new dbdet_curvelet();
//      
//      // decide between ABXY or ABYX using distance
//      if (vgl_distance(eB->pt, eX->pt) < vgl_distance(eB->pt, eY->pt)){
//        quad->push_back(eA, eAcb);
//        quad->push_back(eB, eBcb);
//        quad->push_back(eX, 0);      //no curve bundles for this edgel
//        quad->push_back(eY, 0);      //no curve bundles for this edgel
//      }
//      else {
//        quad->push_back(eA, eAcb);
//        quad->push_back(eB, eBcb);
//        quad->push_back(eY, 0);      //no curve bundles for this edgel
//        quad->push_back(eX, 0);      //no curve bundles for this edgel
//      }
//
//      //record this quad
//      add_a_quadruplet(quad, 0, 1, t1, t2);
//    }
//  }
//
//  // RULE 2. XAB + YAB = XYAB or YXAB 
//  if (t1->edgel_chain[1] == t2->edgel_chain[1] && 
//      t1->edgel_chain[2] == t2->edgel_chain[2])
//  { 
//    dbdet_edgel* eX = t1->edgel_chain[0];
//    dbdet_edgel* eY = t2->edgel_chain[0];
//    dbdet_edgel* eA = t1->edgel_chain[1];//common edgel
//    dbdet_edgel* eB = t1->edgel_chain[2];//common edgel
//    
//    //compute intersection of curve bundles for each perturbation of the common edgels
//    dbdet_curve_bundle *eAcb=0, *eBcb=0;
//    if (t1->curve_bundles[1] && t2->curve_bundles[1])
//      eAcb = new dbdet_curve_bundle(eA, t1->curve_bundles[1], 
//                                        t2->curve_bundles[1] );
//
//    if (t1->curve_bundles[2] && t2->curve_bundles[2])
//      eBcb = new dbdet_curve_bundle(eB, t1->curve_bundles[2], 
//                                        t2->curve_bundles[2] );
//
//    //The participating curve bundles at at least one of the common edgels have to be defined in order to form a quad
//    if (!eAcb && !eBcb)
//      return;
//
//    //if any of the perturbations produced a valid bundle, form an edgel quad
//    if ((eAcb && eAcb->any_valid_bundles()) || 
//        (eBcb && eBcb->any_valid_bundles()))
//    {
//      dbdet_curvelet* quad = new dbdet_curvelet();
//
//      // decide between XYAB or YXAB using distance
//      if (vgl_distance(eA->pt, eX->pt) < vgl_distance(eA->pt, eY->pt)){
//        quad->push_back(eY, 0);      //no curve bundles for this edgel
//        quad->push_back(eX, 0);      //no curve bundles for this edgel
//        quad->push_back(eA, eAcb);
//        quad->push_back(eB, eBcb);
//      }
//      else {
//        quad->push_back(eX, 0);      //no curve bundles for this edgel
//        quad->push_back(eY, 0);      //no curve bundles for this edgel
//        quad->push_back(eA, eAcb);
//        quad->push_back(eB, eBcb);
//      }
//
//      //record this quad
//      add_a_quadruplet(quad, 2, 3, t1, t2);
//    }
//  }
//
//  // RULE 3. XAB + ABY = XABY 
//  if (t1->edgel_chain[1] == t2->edgel_chain[0] && 
//      t1->edgel_chain[2] == t2->edgel_chain[1])
//  { 
//    dbdet_edgel* eX = t1->edgel_chain[0];
//    dbdet_edgel* eA = t1->edgel_chain[1];//common edgel
//    dbdet_edgel* eB = t1->edgel_chain[2];//common edgel
//    dbdet_edgel* eY = t2->edgel_chain[2];
//
//    num_ABC_BCD++;
//
//    //compute intersection of curve bundles for each perturbation of the common edgels
//    dbdet_curve_bundle *eAcb=0, *eBcb=0;
//    if (t1->curve_bundles[1] && t2->curve_bundles[0])
//      eAcb = new dbdet_curve_bundle(eA, t1->curve_bundles[1], 
//                                        t2->curve_bundles[0] );
//
//    if (t1->curve_bundles[2] && t2->curve_bundles[1])
//      eBcb = new dbdet_curve_bundle(eB, t1->curve_bundles[2], 
//                                        t2->curve_bundles[1] );
//
//    //The participating curve bundles at at least one of the common edgels have to be defined in order to form a quad
//    if (!eAcb && !eBcb)
//      return;
//
//    num_ABC_BCD_bundles_ok++;
//
//    //if any of the perturbations produced a valid bundle, form an edgel quad
//    if ((eAcb && eAcb->any_valid_bundles()) || 
//        (eBcb && eBcb->any_valid_bundles()))
//    {
//      dbdet_curvelet* quad = new dbdet_curvelet();
//      quad->push_back(eX, 0);      //no curve bundles for this edgel
//      quad->push_back(eA, eAcb);
//      quad->push_back(eB, eBcb);
//      quad->push_back(eY, 0);      //no curve bundles for this edgel
//
//      //record this triplet
//      add_a_quadruplet(quad, 1, 2, t1, t2);
//    }
//    else
//      num_ABC_BCD_intersection_invalid++;
//  }
//
//  // RULE 4. ABX + AYB = AYBX 
//  if (t1->edgel_chain[0] == t2->edgel_chain[0] && 
//      t1->edgel_chain[1] == t2->edgel_chain[2])
//  { 
//    dbdet_edgel* eA = t1->edgel_chain[0];//common edgel
//    dbdet_edgel* eB = t1->edgel_chain[1];//common edgel
//    dbdet_edgel* eX = t1->edgel_chain[2];
//    dbdet_edgel* eY = t2->edgel_chain[1];
//
//    //compute intersection of curve bundles for each perturbation of the common edgels
//    dbdet_curve_bundle *eAcb=0, *eBcb=0;
//    if (t1->curve_bundles[0] && t2->curve_bundles[0])
//      eAcb = new dbdet_curve_bundle(eA, t1->curve_bundles[0], 
//                                        t2->curve_bundles[0] );
//
//    if (t1->curve_bundles[1] && t2->curve_bundles[2])
//      eBcb = new dbdet_curve_bundle(eB, t1->curve_bundles[1], 
//                                        t2->curve_bundles[2] );
//
//    //The participating curve bundles at at least one of the common edgels have to be defined in order to form a quad
//    if (!eAcb && !eBcb)
//      return;
//
//    //if any of the perturbations produced a valid bundle, form an edgel quad
//    if ((eAcb && eAcb->any_valid_bundles()) || 
//        (eBcb && eBcb->any_valid_bundles()))
//    {
//      dbdet_curvelet* quad = new dbdet_curvelet();
//      quad->push_back(eA, eAcb);
//      quad->push_back(eY, 0);      //no curve bundles for this edgel
//      quad->push_back(eB, eBcb);
//      quad->push_back(eX, 0);      //no curve bundles for this edgel
//
//      //record this quad
//      add_a_quadruplet(quad, 0, 2, t1, t2);
//    }
//  }
//
//  // RULE 5. XAB + AYB = XAYB 
//  if (t1->edgel_chain[1] == t2->edgel_chain[0] && 
//      t1->edgel_chain[2] == t2->edgel_chain[2])
//  { 
//    dbdet_edgel* eX = t1->edgel_chain[0];
//    dbdet_edgel* eA = t1->edgel_chain[1];//common edgel
//    dbdet_edgel* eB = t1->edgel_chain[2];//common edgel
//    dbdet_edgel* eY = t2->edgel_chain[1];
//
//    //compute intersection of curve bundles for each perturbation of the common edgels
//    dbdet_curve_bundle *eAcb=0, *eBcb=0;
//    if (t1->curve_bundles[1] && t2->curve_bundles[0])
//      eAcb = new dbdet_curve_bundle(eA, t1->curve_bundles[1], 
//                                        t2->curve_bundles[0] );
//
//    if (t1->curve_bundles[2] && t2->curve_bundles[2])
//      eBcb = new dbdet_curve_bundle(eB, t1->curve_bundles[2], 
//                                        t2->curve_bundles[2] );
//
//    //The participating curve bundles at at least one of the common edgels have to be defined in order to form a quad
//    if (!eAcb && !eBcb)
//      return;
//
//    //if any of the perturbations produced a valid bundle, form an edgel quad
//    if ((eAcb && eAcb->any_valid_bundles()) || 
//        (eBcb && eBcb->any_valid_bundles()))
//    {
//      dbdet_curvelet* quad = new dbdet_curvelet();
//      quad->push_back(eX, 0);      //no curve bundles for this edgel
//      quad->push_back(eA, eAcb);
//      quad->push_back(eY, 0);      //no curve bundles for this edgel
//      quad->push_back(eB, eBcb);
//
//      //record this quad
//      add_a_quadruplet(quad, 1, 3, t1, t2);
//    }
//  }
//
//  // RULE 6. AXB + AYB = AXYB or AYXB 
//  if (t1->edgel_chain[0] == t2->edgel_chain[0] && 
//      t1->edgel_chain[2] == t2->edgel_chain[2])
//  { 
//    dbdet_edgel* eA = t1->edgel_chain[0];//common edgel
//    dbdet_edgel* eX = t1->edgel_chain[1];
//    dbdet_edgel* eB = t1->edgel_chain[2];//common edgel
//    dbdet_edgel* eY = t2->edgel_chain[1];
//
//    //compute intersection of curve bundles for each perturbation of the common edgels
//    dbdet_curve_bundle *eAcb=0, *eBcb=0;
//    if (t1->curve_bundles[0] && t2->curve_bundles[0])
//      eAcb = new dbdet_curve_bundle(eA, t1->curve_bundles[0], 
//                                        t2->curve_bundles[0] );
//
//    if (t1->curve_bundles[2] && t2->curve_bundles[2])
//      eBcb = new dbdet_curve_bundle(eB, t1->curve_bundles[2], 
//                                        t2->curve_bundles[2] );
//
//    //The participating curve bundles at at least one of the common edgels have to be defined in order to form a quad
//    if (!eAcb && !eBcb)
//      return;
//
//    //if any of the perturbations produced a valid bundle, form an edgel quad
//    if ((eAcb && eAcb->any_valid_bundles()) || 
//        (eBcb && eBcb->any_valid_bundles()))
//    {
//      dbdet_curvelet* quad = new dbdet_curvelet();
//
//      // decide between AXYB or AYXB using distance
//      if (vgl_distance(eA->pt, eX->pt) < vgl_distance(eA->pt, eY->pt)){
//        quad->push_back(eA, eAcb);
//        quad->push_back(eX, 0);      //no curve bundles for this edgel
//        quad->push_back(eY, 0);      //no curve bundles for this edgel
//        quad->push_back(eB, eBcb);
//      }
//      else {
//        quad->push_back(eA, eAcb);
//        quad->push_back(eY, 0);      //no curve bundles for this edgel
//        quad->push_back(eX, 0);      //no curve bundles for this edgel
//        quad->push_back(eB, eBcb);
//      }
//
//      //record this quad
//      add_a_quadruplet(quad, 0, 3, t1, t2);
//    }
//  }
//
//}
//
////: add a newly formed edgel quadruplet to all relevant lists
//void dbdet_sel::add_a_quadruplet(dbdet_curvelet* quad, unsigned pos1, unsigned pos2,
//                                  dbdet_curvelet* t1, dbdet_curvelet* t2)
//{
//  //need to do a quad existence check first because this quad might already 
//  //have formed due to some other rule
//  dbdet_curvelet* t1_prev_quad = t1->exists(quad);
//  dbdet_curvelet* t2_prev_quad = t2->exists(quad);
//
//  //if it exists in t1
//  if (t1_prev_quad)
//  {
//    //if it already exists check to see if the common edgel already has a curve bundle
//    if (!t1_prev_quad->curve_bundles[pos1]){
//      t1_prev_quad->curve_bundles[pos1] = quad->curve_bundles[pos1];
//      quad->curve_bundles[pos1] = 0;
//    }
//    if (!t1_prev_quad->curve_bundles[pos2]){
//      t1_prev_quad->curve_bundles[pos2] = quad->curve_bundles[pos2];
//      quad->curve_bundles[pos2] = 0;
//    }
//    delete quad; //this quad is not required anymore
//
//    //if t2 is not linked to this quad, make a link
//    if (!t2_prev_quad)
//      t2->add_larger_curvelet(t1_prev_quad);
//  }
//  else if (t2_prev_quad) //if it exists in t2
//  {  
//    //if it already exists check to see if the common edgel already has a curve bundle
//    if (!t2_prev_quad->curve_bundles[pos1]){
//      t2_prev_quad->curve_bundles[pos1] = quad->curve_bundles[pos1];
//      quad->curve_bundles[pos1] = 0;
//    }
//    if (!t2_prev_quad->curve_bundles[pos2]){
//      t2_prev_quad->curve_bundles[pos2] = quad->curve_bundles[pos2];
//      quad->curve_bundles[pos2] = 0;
//    }
//    delete quad; //this is not required anymore
//
//    //if t1 is not linked to this quad, make a link
//    if (!t1_prev_quad)
//      t1->add_larger_curvelet(t2_prev_quad);
//  }
//  else { 
//    //brand new quad
//
//    //add links from the edgels to this triplet
//    quad->add_links_from_participating_edgels();
//
//    //also make links from the participating pairs
//    t1->add_larger_curvelet(quad);
//    t2->add_larger_curvelet(quad);
//
//    //and finally, add it to the edge quads list
//    quads.push_back(quad);
//  }
//}

////: form an edgel grouping by adding a pair to the existing curvelet
//template <class curve_model>
//void dbdet_sel<curve_model>::form_an_edgel_grouping(dbdet_edgel* eA, dbdet_curvelet* cvlet, dbdet_curvelet* pair)
//{
//  //1) first determine the index of the current edgel in the groupings
//  //   Also determine if this pair is already contained in the current cvlet.
//
//  unsigned cvlet_ind = -1; //index of the ref edgel in the cvlet
//  unsigned pair_ind = -1;  //index of the ref edgel in the pair
//  dbdet_edgel* eN = 0; //new edgel to be added
//  
//  if (pair->edgel_chain[0]==eA) {pair_ind = 0; eN = pair->edgel_chain[1];}
//  else                          {pair_ind = 1; eN = pair->edgel_chain[0];}                    
//
//  for (unsigned i=0; i<cvlet->edgel_chain.size(); i++){
//    if (cvlet->edgel_chain[i]==eA)
//      cvlet_ind = i;
//    else if (cvlet->edgel_chain[i]==eN){
//      cvlet_ind = -1; //duplicate edgel pair
//      break;
//    }
//  }
//  
//  if (cvlet_ind==-1) // this pair is already contained in the current curvelet
//    return;
//
//  //2) if not, compute intersection of curve bundles for each perturbation of the common edgel
//  dbdet_curve_bundle *eAcb=0;
//  if (cvlet->curve_bundles[cvlet_ind] && pair->curve_bundles[pair_ind])
//    eAcb = new dbdet_curve_bundle(eA, cvlet->curve_bundles[cvlet_ind], 
//                                      pair->curve_bundles[pair_ind]);
//
//  //if any of the perturbations produced a valid bundle, form an new edgel grouping
//  if (eAcb && eAcb->any_valid_bundles())
//  {
//    dbdet_curvelet* new_cvlet = new dbdet_curvelet();
//
//    //decide the placement of the new edgel in the cvlet by computing distances
//    unsigned new_edge_ind = -1;
//    unsigned new_ref_edge_ind = -1;
//    if (pair_ind==1){ //new edgel goes before the ref edgel
//      for (unsigned i=0; i<cvlet_ind; i++){
//        if (vgl_distance(eA->pt, eN->pt) > vgl_distance(eA->pt, cvlet->edgel_chain[i]->pt)){
//          new_edge_ind = i;
//          break;
//        }
//      }
//      if (new_edge_ind==-1) //new edgel gets added in front of the ref edgel
//        new_edge_ind = cvlet_ind;
//
//      new_ref_edge_ind = cvlet_ind+1; //new index of the reference edgel
//    }
//    else { //new edgel goes after the ref edgel
//      for (unsigned i=cvlet_ind+1; i<cvlet->edgel_chain.size(); i++){
//        if (vgl_distance(eA->pt, eN->pt) < vgl_distance(eA->pt, cvlet->edgel_chain[i]->pt)){
//          new_edge_ind = i;
//          break;
//        }
//      }
//      if (new_edge_ind==-1) //new edgel gets added at the end
//        new_edge_ind = cvlet->edgel_chain.size();
//
//      new_ref_edge_ind = cvlet_ind; //new index of the reference edgel
//    }
//
//    //construct the ordered edgel chain for the new curvelet
//    for (unsigned i=0; i<cvlet->edgel_chain.size(); i++){
//      if (i==new_edge_ind) //add the new edgel in the middle
//        new_cvlet->push_back(eN, 0);
//
//      if (i==cvlet_ind)
//        new_cvlet->push_back(eA, eAcb); //new curve bundle
//      else
//        new_cvlet->push_back(cvlet->edgel_chain[i], 0); //no curve bundles for this edgel
//    }
//    if (new_edge_ind==cvlet->edgel_chain.size()) //new edgels that go to the end
//      new_cvlet->push_back(eN, 0);
//
//    //record this curvelet
//    add_a_curvelet(new_cvlet, new_ref_edge_ind, cvlet);
//  }
//}

////: add a newly formed curvelet to all the relevant lists
//template <class curve_model>
//void dbdet_sel<curve_model>::add_a_curvelet(dbdet_curvelet* new_cvlet, unsigned ref_edge_ind, dbdet_curvelet* old_cvlet)
//{
//  dbdet_edgel* ref_edgel = new_cvlet->edgel_chain[ref_edge_ind];
//  unsigned order = new_cvlet->edgel_chain.size();
//
//  //need to do an existence check first because this curvelet might already 
//  //have formed due to some other rule
//  dbdet_curvelet* prev_cvlet = 0;
//
//  //go over all groupings of this order formed by the ref edgel
//  if (ref_edgel->local_curvelets.size()>order-2)
//  {
//    curvelet_list_iter prev_it = ref_edgel->local_curvelets[order-2].begin();
//    for (; prev_it != ref_edgel->local_curvelets[order-2].end(); prev_it++){
//      dbdet_curvelet* lcv = (*prev_it);
//
//      bool found = true;
//      for (unsigned j=0; j<lcv->edgel_chain.size(); j++){
//        if (lcv->edgel_chain[j] != new_cvlet->edgel_chain[j]){
//          found = false;
//          break;
//        }
//      }
//
//      if (found){
//        prev_cvlet = lcv;
//        break;
//      }
//    }
//  }
//
//  //if it exists
//  if (prev_cvlet)
//  {
//    //if the curvelet already exists check to see if the common edgel already has a curve bundle
//    if (!prev_cvlet->curve_bundles[ref_edge_ind]){
//      prev_cvlet->curve_bundles[ref_edge_ind] = new_cvlet->curve_bundles[ref_edge_ind];
//      new_cvlet->curve_bundles[ref_edge_ind] = 0;
//    }
//    delete new_cvlet; //this curvelet is not required anymore
//  }
//  else { 
//    //brand new curvelet
//
//    //add links from the ref edgel to this curvelet
//    ref_edgel->add_curvelet(new_cvlet);
//
//    //also make links from the participating curvelets
//    //old_cvlet->add_larger_curvelet(new_cvlet);
//
//    //and finally, add it to the edge groupings list
//    //quads.push_back(quad);
//  }
//}

//: form curvelets around the given edgel in a greedy fashion
template <class curve_model>
void dbdet_sel<curve_model>::build_curvelets_greedy_for_edge(dbdet_edgel* eA, unsigned max_size_to_group, bool use_flag, 
                                bool forward, bool centered, bool leading)
{
  // 1) construct a structure to temporarily hold the pairwise-hypotheses
  vcl_vector<sel_hyp> eA_hyps;

  //get the grid coordinates of this edgel
  unsigned const ii = dbdet_round(eA->pt.x());
  unsigned const jj = dbdet_round(eA->pt.y());

  unsigned const rad_sqr = rad_*rad_;
  
  // 2) iterate over the neighboring cells around this edgel
  for (int xx=(int)ii-(int)nrad_; xx<=(int)(ii+nrad_) ; xx++){
    for (int yy=(int)jj-(int)nrad_; yy<=(int)(jj+nrad_) ; yy++){

      if (xx<0 || xx>=(int)ncols_ || yy<0 || yy>=(int)nrows_)
        continue;

      //for all the edgels in its neighborhood
      for (unsigned k=0; k<edgemap_->cell(xx, yy).size(); k++){
        dbdet_edgel* eB = edgemap_->cell(xx, yy)[k];
        if (eB == eA) continue;

        // 3) do a better check of circular radius because the bucketing neighborhood is very coarse
        if ( sqr_length(eA->pt - eB->pt) > rad_sqr) continue;

        // 4) form pair-wise hypotheses
        bool ref_first;
        curve_model* cm = form_a_hypothesis(eA, eB, ref_first, forward, centered, leading);

        if (cm){
          if (cm->bundle_is_valid()){ //if legal, record the hypothesis
            sel_hyp cur_hyp;

            cur_hyp.eN = eB;
            cur_hyp.cm = cm;
            cur_hyp.ref_first = ref_first;
            cur_hyp.d = vgl_distance(eA->pt, eB->pt);
            cur_hyp.flag = false;

            eA_hyps.push_back(cur_hyp);
          }
          else
            delete cm;
        }
      }
    }
  }

  // 5) first sort the pair-wise hyps by distance
  vcl_sort(eA_hyps.begin(), eA_hyps.end(), comp_dist_hyps_less);

  // 6) for each pair-wise hyps formed by this edgel
  for (unsigned h1=0; h1<eA_hyps.size(); h1++)
  {
    //this edgel has been used up already (only if we are using flags)
    if (use_flag && eA_hyps[h1].flag) continue;

    //for the hybrid algorithm (to avoid edgel skips within contour)
    if (use_hybrid_){ if (cId_[eA->id]==cId_[eA_hyps[h1].eN->id]) continue;}

    // 7) initialize a new edgel chain that will grow in a greedy depth-first fashion
    vcl_deque<dbdet_edgel*> cur_edgel_chain; //chain can grow either way

    //insert ref edgel first
    cur_edgel_chain.push_back(eA); 
      
    // 8)initialize the shrinking curve bundle from the cur hypothesis
    curve_model* cur_cm = static_cast<curve_model *>(eA_hyps[h1].cm);

    // 9) attempt to integrate the other pair-wise hypotheses (DEPTH-FIRST SEARCH)
    for (unsigned h2=0; h2<eA_hyps.size(); h2++)
    {
      if (h2==h1){
        //add the second edgel of the reference pair (no need for any computation)
        if (eA_hyps[h1].ref_first) cur_edgel_chain.push_back(eA_hyps[h1].eN);
        else                       cur_edgel_chain.push_front(eA_hyps[h1].eN);

        // check the size of the grouping
        if (cur_edgel_chain.size() >= max_size_to_group)
          break;
        else
          continue;
      }

      // 10) for the others we need to check for consistency with the current hypothesis
      //     (i.e., compute intersection of curve bundles)
      curve_model* new_cm = cur_cm->intersect(static_cast<curve_model *>(eA_hyps[h2].cm));

      // 11) if the intersection is valid, we can add this edgel to the grouping
      if (new_cm->bundle_is_valid())
      {
        //reassign the curve bundle for the growing grouping
        if (cur_cm != eA_hyps[h1].cm)
          delete cur_cm;
        cur_cm = new_cm;

        // 12) add the new edgel to the growing edgel chain
        if (eA_hyps[h2].ref_first) cur_edgel_chain.push_back(eA_hyps[h2].eN);
        else                       cur_edgel_chain.push_front(eA_hyps[h2].eN);

        //flag this hypothesis
        eA_hyps[h2].flag = true;
      }
      else
        delete new_cm; //delete this cb because it is not needed

      // 13) check the size of the grouping
      if (cur_edgel_chain.size() >= max_size_to_group)
        break;

    }

    // 14) form a new curvelet and assign the curve bundle and the edgel chain to it...
    //     ...if it passes a few tests : 
    //                             (a) if it doesn't already exist and 
    //                             (b) its a reasonable fit and
    //                             (c) its relatively symmetric
    if (cur_edgel_chain.size()>2 && 
        !curvelet_map_.does_curvelet_exist(eA, cur_edgel_chain) && 
        cur_cm->curve_fit_is_reasonable(cur_edgel_chain, eA, dpos_))// &&
        //curvelet_is_balanced(eA, cur_edgel_chain))
    {
      dbdet_curvelet* new_cvlet = new dbdet_curvelet(eA, cur_cm, cur_edgel_chain, forward);

      //compute curvelet quality
      new_cvlet->compute_properties(rad_, token_len_);

      //add links from the ref edgel to this curvelet
      curvelet_map_.add_curvelet(new_cvlet);
    }
    else {
      //delete cur_cm if no curvelet is formed
      if (cur_cm != eA_hyps[h1].cm)
        delete cur_cm;
    }
  }

  // 15) now delete all the hyps for this edgel, before moving on to a new edgel
  for (unsigned h1=0; h1<eA_hyps.size(); h1++)
    delete eA_hyps[h1].cm;

  eA_hyps.clear();
}

//: form an edgel grouping from an ordered list of edgels
template <class curve_model>
dbdet_curvelet* 
dbdet_sel<curve_model>::form_an_edgel_grouping(dbdet_edgel* ref_e, 
    vcl_deque<dbdet_edgel*> &edgel_chain, 
    bool forward,  bool centered, bool leading)
{
  //1) Go over the edgels in the chain and attempt to form a curvelet from it
  curve_model* chain_cm=0;

  // 2) form an ordered list based on distance from ref_edgel
  vcl_map<double, unsigned > dist_order; // ordered list of edgels (closest to furthest)
  for (unsigned i=0; i< edgel_chain.size(); i++)
  {
    if (edgel_chain[i] == ref_e)   continue;
    dist_order.insert(vcl_pair<double, unsigned>(vgl_distance(ref_e->pt, edgel_chain[i]->pt), i));
  } 

  vcl_map<double, unsigned >::iterator it = dist_order.begin();
  for (; it!=dist_order.end(); it++)
  {
    dbdet_edgel* cur_e = edgel_chain[it->second];

    //form a pairwise curve bundle between the current edgel and the ref edgel
    bool ref_first; //this is not used here
    curve_model* cm = form_a_hypothesis(ref_e, cur_e, ref_first, forward, centered, leading);
    
    //if the bundle is legal intersect it with the current group's bundle
    if (cm && cm->bundle_is_valid())
    {
      //if this is the first pair, record it as the bundle of the grouping
      if (!chain_cm)
        chain_cm = cm;
      else { 
        //intersect it with the existing bundle
        curve_model* new_cm = chain_cm->intersect(cm);
        delete cm; //no use for this anymore
        
        //if this intersection is valid, make it the chain's curve bundle
        if (new_cm->bundle_is_valid()){
          delete chain_cm;   //no use for this since it will be replaced
          chain_cm = new_cm; //replace it
        }
        else {
          //delete any existing curve bundles
          delete new_cm;

          if (chain_cm)  delete chain_cm;
          chain_cm = 0;

          break; //break the for loop (no curve bundle possible)
        }
      }
    }
    else {
      // the entire curvelet is invalid if anyone pairwise curve bundle is invalid
      //delete any existing curve bundles
      if (cm)        delete cm;
      if (chain_cm)  delete chain_cm;
      chain_cm = 0;

      break; //break the for loop (no curve bundle possible)
    }
  }

  dbdet_curvelet* new_cvlet=0;
  //if all the pairwise bundles intersect,  
  //form a new curvelet and assign the curve bundle and the edgel chain to it...
  if (chain_cm)
  {
    chain_cm->curve_fit_is_reasonable(edgel_chain, ref_e, dpos_); //this is just to find the centroid
  
    new_cvlet = new dbdet_curvelet(ref_e, chain_cm, forward);

    //add edgels
    for (unsigned i=0; i<edgel_chain.size(); i++)
      new_cvlet->push_back(edgel_chain[i]);

    //compute curvelet quality
    new_cvlet->compute_properties(rad_, token_len_);

    //Note:: since this curvelet could be temporary, don't add it to the edgel's list of curvelets
    //curvelet_map_.add_curvelet(new_cvlet);
  }

  return new_cvlet;
}
