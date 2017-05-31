// This is brcv/dbdet/algo/dbdet_postprocess_contours.cxx

//:
// \file

#include "dbdet_postprocess_contours.h"

#include "dbdet_EMD.h"
#include "bgld_arc_algo.h"
#include "bgld_curve_smoothing.h"
#include <mbl/mbl_stats_1d.h>
#include <vcl_algorithm.h>

//: post process to break contours at high curvature points
void post_process_based_on_curvature(dbdet_curve_fragment_graph& curve_frag_graph, double k_thresh)
{
  //since the chains are being fragmented we cannot use the iterators to 
  //find the end of the original fragment list
  unsigned orig_frag_number = curve_frag_graph.frags.size();

  unsigned cnt=0;
  dbdet_edgel_chain_list_iter f_it = curve_frag_graph.frags.begin();
  for (; f_it != curve_frag_graph.frags.end(); f_it++)
  {
    dbdet_edgel_chain* chain = (*f_it);

    if (++cnt>orig_frag_number) //terminate after we have gone over the original list of fragments
      break;

    vcl_vector<double> ks(chain->edgels.size());
    vcl_vector<unsigned> break_pts;

    for (unsigned j=0; j<chain->edgels.size(); j++)
    {
      if (chain->edgels.size()<3){
        ks[j] = 0.0;
      }
      else{
        unsigned k1; //index of the first of 3 pts
        if (j<1)                            k1=0;
        else if  (j>chain->edgels.size()-2) k1=chain->edgels.size()-3;
        else                                k1 = j-1;

        //compute curvature
        ks[j] = 1/bgld_arc_algo::compute_arc_radius_from_three_points(chain->edgels[k1]->pt, chain->edgels[k1+1]->pt, chain->edgels[k1+2]->pt);
        
        if (vcl_fabs(ks[j])>k_thresh) 
          break_pts.push_back(j); //needs to be broken here
      }
    }

    //break up this contour
    if (break_pts.size()>0)
    {
      unsigned last_j = 0;

      for (unsigned i=0; i<break_pts.size(); i++){
        //construct new chains from the sub fragments if they are longer than 3

        if ((break_pts[i]-last_j)>3){
          dbdet_edgel_chain* new_chain = new dbdet_edgel_chain();

          for (unsigned j=last_j; j<break_pts[i]-1; j++)
            new_chain->push_back(chain->edgels[j]);

          //add the new fragment
          curve_frag_graph.insert_fragment(new_chain);

          last_j = break_pts[i];
        }
      }

      //last segment
      if ((chain->edgels.size()-last_j)>3)
      {
        dbdet_edgel_chain* new_chain = new dbdet_edgel_chain();

        for (unsigned j=last_j; j<chain->edgels.size(); j++)
          new_chain->push_back(chain->edgels[j]);

        //add the new fragment
        curve_frag_graph.insert_fragment(new_chain);
      }

      //remove the original fragment
      if (f_it==curve_frag_graph.frags.begin()){
        curve_frag_graph.remove_fragment(chain);
        f_it = curve_frag_graph.frags.begin();
      }
      else {
        f_it--;
        curve_frag_graph.remove_fragment(chain);
      }
    }
  }
}

void appearance_based_post_processing(dbdet_curve_fragment_graph& curve_frag_graph, unsigned win_len, double adap_thresh)
{
  //since the chains are being fragmented we cannot use the iterators to 
  //find the end of the original fragment list
  unsigned orig_frag_number = curve_frag_graph.frags.size();

  unsigned cnt=0;
  dbdet_edgel_chain_list_iter f_it = curve_frag_graph.frags.begin();
  for (; f_it != curve_frag_graph.frags.end(); f_it++)
  {
    dbdet_edgel_chain* chain = (*f_it);

    if (++cnt>orig_frag_number) //terminate after we have gone over the original list of fragments
      break;

    vcl_deque<double> Lwin;
    vcl_deque<double> Rwin;

    vcl_vector<unsigned> break_pts;

    for (unsigned j=0; j<chain->edgels.size(); j++)
    {
      //add value to the buffer
      Lwin.push_back(chain->edgels[j]->left_app->value());
      Rwin.push_back(chain->edgels[j]->right_app->value());

      if (Lwin.size()>win_len) Lwin.pop_front();
      if (Rwin.size()>win_len) Rwin.pop_front();
      
      mbl_stats_1d Ldata, Rdata;

      for (unsigned i=0; i< vcl_min((size_t)win_len, Lwin.size()); i++) Ldata.obs(Lwin[i]);
      for (unsigned i=0; i< vcl_min((size_t)win_len, Rwin.size()); i++) Rdata.obs(Rwin[i]);

      if (vcl_fabs(Ldata.mean()-Rdata.mean())< adap_thresh*(Ldata.sd()+Rdata.sd())) 
      {
        break_pts.push_back(j); //needs to be broken here
        Lwin.clear(); //reset buffer
        Rwin.clear();
      }
    }

    //break up this contour
    if (break_pts.size()>0)
    {

      unsigned last_j = 0;

      for (unsigned i=0; i<break_pts.size(); i++){
        //construct new chains from the sub fragments if they are longer than 3

        if ((break_pts[i]-last_j)>3){
          dbdet_edgel_chain* new_chain = new dbdet_edgel_chain();

          for (unsigned j=last_j; j<break_pts[i]-1; j++)
            new_chain->push_back(chain->edgels[j]);

          //add the new fragment
          curve_frag_graph.insert_fragment(new_chain);

          last_j = break_pts[i];
        }
      }

      //last segment
      if ((chain->edgels.size()-last_j)>3)
      {
        dbdet_edgel_chain* new_chain = new dbdet_edgel_chain();

        for (unsigned j=last_j; j<chain->edgels.size(); j++)
          new_chain->push_back(chain->edgels[j]);

        //add the new fragment
        curve_frag_graph.insert_fragment(new_chain);
      }

      //remove the original fragment
      if (f_it==curve_frag_graph.frags.begin()){
        curve_frag_graph.remove_fragment(chain);
        f_it = curve_frag_graph.frags.begin();
      }
      else {
        f_it--;
        curve_frag_graph.remove_fragment(chain);
      }
    }
  }

}

double avg_strength(dbdet_edgel_chain* chain)
{
  //collect other statistics for the contours
  mbl_stats_1d data;

  // 1) compute the average edge strength of an edgel chain
  for (unsigned j=0; j<chain->edgels.size(); j++)
    data.obs(chain->edgels[j]->strength);

  return data.mean();
}

double mean_contrast(dbdet_edgel_chain* chain)
{
  if (dynamic_cast<dbdet_intensity*>(chain->edgels[0]->left_app))
  {
    mbl_stats_1d Ldata, Rdata;

    // compute mean   
    Ldata.clear();  Rdata.clear();
    for (unsigned k=0; k<chain->edgels.size(); k++){
      Ldata.obs(chain->edgels[k]->left_app->value());
      Rdata.obs(chain->edgels[k]->right_app->value());
    }

    return vcl_fabs(Ldata.mean() - Rdata.mean());
  }
  else if (dynamic_cast<dbdet_color*>(chain->edgels[0]->left_app))
  {
    mbl_stats_1d LLdata, Ladata, Lbdata;
    mbl_stats_1d RLdata, Radata, Rbdata;

    // compute mean   
    LLdata.clear(); Ladata.clear(); Lbdata.clear();
    RLdata.clear(); Radata.clear(); Rbdata.clear();

    for (unsigned k=0; k<chain->edgels.size(); k++){
      LLdata.obs(((dbdet_color*)chain->edgels[k]->left_app)->c1);
      Ladata.obs(((dbdet_color*)chain->edgels[k]->left_app)->c2);
      Lbdata.obs(((dbdet_color*)chain->edgels[k]->left_app)->c3);
      
      RLdata.obs(((dbdet_color*)chain->edgels[k]->right_app)->c1);
      Radata.obs(((dbdet_color*)chain->edgels[k]->right_app)->c2);
      Rbdata.obs(((dbdet_color*)chain->edgels[k]->right_app)->c3);
    }

    double dL = LLdata.mean()-RLdata.mean();
    double da = Ladata.mean()-Radata.mean();
    double db = Lbdata.mean()-Rbdata.mean();

    return vcl_sqrt(dL*dL+da*da+db*db);
  }
  else if (dynamic_cast<dbdet_gray_signature*>(chain->edgels[0]->left_app))
  {
    dbdet_signature Lcomb, Rcomb;

    //compute aggregate histograms
    for (unsigned k=0; k<chain->edgels.size(); k++){
      Lcomb += ((dbdet_gray_signature*)chain->edgels[k]->left_app)->sig;
      Rcomb += ((dbdet_gray_signature*)chain->edgels[k]->right_app)->sig;
    }

    return Lcomb-Rcomb;
  }
  else 
    return 0.0;

}

double left_app_std(dbdet_edgel_chain* chain)
{
  if (dynamic_cast<dbdet_intensity*>(chain->edgels[0]->left_app))
  {
    mbl_stats_1d data;

    // 1) compute the average edge strength of an edgel chain
    for (unsigned j=0; j<chain->edgels.size(); j++)
      data.obs(chain->edgels[j]->left_app->value());

    return data.sd();
  }
  else 
    return 0.0; //don't know what to do yet
}

double right_app_std(dbdet_edgel_chain* chain)
{
  if (dynamic_cast<dbdet_intensity*>(chain->edgels[0]->left_app))
  {
    mbl_stats_1d data;

    // 1) compute the average edge strength of an edgel chain
    for (unsigned j=0; j<chain->edgels.size(); j++)
      data.obs(chain->edgels[j]->right_app->value());

    return data.sd();
  }
  else
    return 0.0; //don't know what to do yet
}

double avg_d2f(dbdet_edgel_chain* chain)
{
  //collect other statistics for the contours
  mbl_stats_1d data;

  // 1) compute the average edge strength of an edgel chain
  for (unsigned j=0; j<chain->edgels.size(); j++)
    data.obs(vcl_fabs(chain->edgels[j]->deriv));

  return data.mean();
}

void compute_curvatures(dbdet_edgel_chain* chain, mbl_stats_1d &ks)
{
  if (chain->edgels.size()<3)
  {
    ks.obs(0.0); ks.obs(0.0);
    return;
  }

  //create a polyline out of the edgel chain
  vcl_vector<vgl_point_2d<double> > pts;
  pts.reserve(chain->edgels.size());
  for (unsigned j=0; j<chain->edgels.size(); j++)
    pts.push_back(chain->edgels[j]->pt);

  // smooth this contour
  bgld_csm(pts, 1.0, 1);

  //compute curvatures
  for (unsigned j=0; j<pts.size(); j++)
  {
    unsigned k1; //index of the first of 3 pts
    if (j<1)                  k1=0;
    else if  (j>pts.size()-2) k1=pts.size()-3;
    else                      k1 = j-1;

    //compute curvature
    ks.obs(vcl_fabs(1/bgld_arc_algo::compute_arc_radius_from_three_points(pts[k1], pts[k1+1], pts[k1+2])));
  }
}

double avg_curvature(dbdet_edgel_chain* chain)
{
  mbl_stats_1d ks;
  compute_curvatures(chain, ks);

  // compute the average curvature
  return ks.mean();
}

double max_curvature(dbdet_edgel_chain* chain)
{
  mbl_stats_1d ks;
  compute_curvatures(chain, ks);

  // compute the average curvature
  return ks.max();
}

//: Prune contours on the basis of the length of the curve
void prune_contours_by_length(dbdet_curve_fragment_graph& curve_frag_graph, 
                              double len_thresh)
{
  //TODO: 
  //  if a prune results in a simple connectivity between the fragments,
  //  they need to be spliced

  vcl_vector<dbdet_edgel_chain*> chains_to_del;

  dbdet_edgel_chain_list_iter f_it = curve_frag_graph.frags.begin();
  for (; f_it != curve_frag_graph.frags.end(); f_it++)
  {
    dbdet_edgel_chain* chain = (*f_it);

    //A) Apply the length threshold
    if (chain->edgels.size()<len_thresh){
      chains_to_del.push_back(chain);
      continue;
    }
  }

  //now delete the contours from the list
  for (unsigned i=0; i<chains_to_del.size(); i++)
    curve_frag_graph.remove_fragment(chains_to_del[i]);

}

//: Prune contours on the basis of the sum strength of the curve
void prune_contours_by_strength(dbdet_curve_fragment_graph& curve_frag_graph, 
                              double strength_thresh)
{
  //TODO: 
  //  if a prune results in a simple connectivity between the fragments,
  //  they need to be spliced

  vcl_vector<dbdet_edgel_chain*> chains_to_del;

  dbdet_edgel_chain_list_iter f_it = curve_frag_graph.frags.begin();
  for (; f_it != curve_frag_graph.frags.end(); f_it++)
  {
    dbdet_edgel_chain* chain = (*f_it);

    //collect other statistics for the contours
    mbl_stats_1d Ldata, Rdata;

    // 1) compute the average edge strength of an edgel chain
    for (unsigned j=0; j<chain->edgels.size(); j++)
      Ldata.obs(chain->edgels[j]->strength);

    double sum_strength = Ldata.sum();

    // B) Apply avg. strength threshold
    if (sum_strength < strength_thresh){
      chains_to_del.push_back(chain);
      continue;
    }
  }

  //now delete the contours from the list
  for (unsigned i=0; i<chains_to_del.size(); i++)
    curve_frag_graph.remove_fragment(chains_to_del[i]);

}
//: Prune contours on the basis of average contrast and length of the curve
void prune_contours(dbdet_curve_fragment_graph& curve_frag_graph, 
                    double len_thresh, 
                    double strength_thresh, 
                    double contrast_thresh, 
                    double adap_thresh_fac, 
                    double d2f_thresh, 
                    double k_thresh)
{
  //TODO: 
  //  if a prune results in a simple connectivity between the fragments,
  //  they need to be spliced

  vcl_vector<dbdet_edgel_chain*> chains_to_del;

  dbdet_edgel_chain_list_iter f_it = curve_frag_graph.frags.begin();
  for (; f_it != curve_frag_graph.frags.end(); f_it++)
  {
    dbdet_edgel_chain* chain = (*f_it);

    //A) Apply the length threshold
    if (chain->edgels.size()<len_thresh){
      chains_to_del.push_back(chain);
      continue;
    }

    //collect other statistics for the contours
    mbl_stats_1d Ldata, Rdata;

    // 1) compute the average edge strength of an edgel chain
    for (unsigned j=0; j<chain->edgels.size(); j++)
      Ldata.obs(chain->edgels[j]->strength);

    double avg_strength = Ldata.mean();

    // B) Apply avg. strength threshold
    if (avg_strength<strength_thresh){
      chains_to_del.push_back(chain);
      continue;
    }

    // 2) compute mean and std intensity of the left and right sides of the contour  
    Ldata.clear();  Rdata.clear();
    for (unsigned k=0; k<chain->edgels.size(); k++){
      Ldata.obs(chain->edgels[k]->left_app->value());
      Rdata.obs(chain->edgels[k]->right_app->value());
    }

    double Lmean = Ldata.mean(); double Rmean = Rdata.mean();
    double Lstd = Ldata.sd(); double Rstd = Rdata.sd(); 

    // C) Apply mean contrast threshold
    if (vcl_fabs(Lmean-Rmean)<contrast_thresh){
      chains_to_del.push_back(chain);
      continue;
    }

    // D) apply adaptive contrast threshold (saliency test)
    if (vcl_fabs(Lmean-Rmean)<adap_thresh_fac*(Lstd+Rstd)){
      chains_to_del.push_back(chain);
      continue;
    }

    // 3) compute average d2f for the contour fragment
    Ldata.clear();  Rdata.clear();
    for (unsigned j=0; j<chain->edgels.size(); j++)
      Ldata.obs(chain->edgels[j]->deriv);

    double d2f_mean = Ldata.mean();

    // E) apply peakiness threshold
    if (vcl_fabs(d2f_mean)<d2f_thresh){
      chains_to_del.push_back(chain);
      continue;
    }

    //-------------------------------------------------------------------
    //create a polyline out of the edgel chain
    vcl_vector<vgl_point_2d<double> > pts;
    pts.reserve(chain->edgels.size());
    for (unsigned j=0; j<chain->edgels.size(); j++)
      pts.push_back(chain->edgels[j]->pt);

    // smooth this contour
    bgld_csm(pts, 1.0, 1);

    vcl_vector<double> ks;
    ks.resize(pts.size());

    //compute curvatures
    for (unsigned j=0; j<pts.size(); j++)
    {
      unsigned k1; //index of the first of 3 pts
      if (j<1)                  k1=0;
      else if  (j>pts.size()-2) k1=pts.size()-3;
      else                      k1 = j-1;

      //compute curvature
      ks[j] = 1/bgld_arc_algo::compute_arc_radius_from_three_points(pts[k1], pts[k1+1], pts[k1+2]);
    }

    // F) Apply curvature threshold
    Ldata.clear();
    for (unsigned j=0; j<ks.size(); j++)
      Ldata.obs(vcl_fabs(ks[j]));

    if (Ldata.mean()>k_thresh){
      chains_to_del.push_back(chain);
      continue;
    }
  }

  //now delete the contours from the list
  for (unsigned i=0; i<chains_to_del.size(); i++)
    curve_frag_graph.remove_fragment(chains_to_del[i]);


}

//: Post process to link distinct but connected curve fragments
void post_process_to_link_fragments(dbdet_curve_fragment_graph& curve_frag_graph)
{  
  // Splice contours if adjacent fragments are unambiguously linked
  for (unsigned i=0; i<curve_frag_graph.size(); i++)
  {
    // if an edgel has exactly one curve fragment before it and one after it
    // merge them
    if (curve_frag_graph.pFrags[i].size()==1 &&
        curve_frag_graph.cFrags[i].size()==1)
    {
      //merge the fragments
      dbdet_edgel_chain* pChain = curve_frag_graph.pFrags[i].front();
      dbdet_edgel_chain* cChain = curve_frag_graph.cFrags[i].front();

      if (pChain==cChain)
        continue;

      pChain->edgels.insert(pChain->edgels.end(), cChain->edgels.begin(), cChain->edgels.end());

      //now remove the second chain and update the curve fragment graph
      curve_frag_graph.pFrags[i].clear();
      curve_frag_graph.pFrags[cChain->edgels.back()->id].clear();
      curve_frag_graph.pFrags[cChain->edgels.back()->id].push_back(pChain);

      curve_frag_graph.remove_fragment(cChain);
    }
  }

}

//: Perform post-processing to remove ambiuous segments between the ends of contour fragments
void post_process_to_remove_ambiguity(dbdet_curve_fragment_graph& /*curve_frag_graph*/)
{

  // Algorithm 2: From the ends of existing curve fragments, track all the paths until we hit another curve segment
  //              Choose the best path between these end points as the correct contour and splice it in.


}


