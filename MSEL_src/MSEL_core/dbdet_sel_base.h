// This is brcv/seg/dbdet/algo/dbdet_sel_base.h
#ifndef dbdet_sel_base_h
#define dbdet_sel_base_h
//:
//\file
//\brief Symbolic edge linking algorithm base
//\author Amir Tamrakar
//\date 09/11/06
//
//\verbatim
//  Modifications
//  Amir Tamrakar 09/05/06   Removed all the other classes to their own files
//  Amir Tamrakar 09/07/06   Templatized this class to enable it to operate on
//                           different curve models
//  Amir Tamrakar 09/11/06   Made a base class that is not templatized so that 
//                           a common smart pointer can be used for all subclasses
//  Amir Tamrakar 09/16/06   Now using the edgemap class instead of edge list so 
//                           there is no need to form edgel neighborhoods
//  Amir Tamrakar 07/20/07   Added upper bound specification for curve bundles
//
//  Amir Tamrakar            Added data structures and algorithms for the DHT algorithm for linking
//
//  Amir Tamrakar            Now passing in cvlet params class instead of all the individual parameters
//                           for curvelet formation
//
//  Amir Tamrakar            Added some post-processing functions to prune contours
//
//  Amir Tamrakar            Curvelets can be formed in 4 ways now: regular anchor centered/ anchor centered bidirectional /
//                           Anchor leading bidirectional / ENO style (anchor leading or trailing but in the same direction)
//
//\endverbatim

#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_array_2d.h>
#include <vcl_vector.h>
#include <vcl_list.h>
#include <vcl_set.h>
#include <vcl_map.h>
#include <vcl_queue.h>
#include <vcl_utility.h>

#include "dbdet_edgel.h"
#include "dbdet_edgemap_sptr.h"
#include "dbdet_edgemap.h"
#include "dbdet_curve_model.h"
#include "dbdet_curvelet.h"
#include "dbdet_curvelet_map.h"
#include "dbdet_edgel_link_graph.h"
#include "dbdet_curve_fragment_graph.h"
#include "dbdet_sel_utils.h"

#include "dbdet_hyp_tree_graph.h"

#include "dbdet_EHT.h"
#include "dbdet_CFTG.h"
//#include "dbdet_sel_knot.h"

//: Base class for the "Symbolic edge linking algorithm" (as described in POCV 06).
//  It contains the following things:
//      a) The algorithm to locally group them into curvelets given
//          -- an estimate of the uncertainty of the edge measurements and
//          -- the regularizing curve model
//
//      b) The algorithm to convert local groups into a link graph
//
//      c) Algorithms to convert a link graph into an image contour topology
//
//  The linking algorithm needs to be able to work with various curve models. 
//  This class is therefore subclasses in a templatized fashion to enable it to work 
//  with various curve models.
//
class dbdet_sel_base : public vbl_ref_count
{
public:

  //: constructor
  dbdet_sel_base(dbdet_edgemap_sptr edgemap, 
                 dbdet_curvelet_map& cvlet_map, 
                 dbdet_edgel_link_graph& edge_link_graph, 
                 dbdet_curve_fragment_graph& curve_frag_graph,
                 dbdet_curvelet_params cvlet_params=dbdet_curvelet_params());  //various parameters

  virtual ~dbdet_sel_base();

  // Access parameters
  unsigned nrows() { return nrows_; }
  unsigned ncols() { return ncols_; }

  double rad() { return rad_; }
  double dtheta() { return dtheta_; }
  double dpos() { return dpos_; }

  unsigned maxN() { return maxN_; }

  //: return a reference to the edgel buckets
  vbl_array_2d<vcl_vector<dbdet_edgel*> > & cells() { return edgemap_->edge_cells; }

  //: return the cvlet map
  dbdet_curvelet_map& cvlet_map() { return curvelet_map_; }

  //: return edge map
  dbdet_edgemap_sptr get_edge_map() { return edgemap_; }

  //: return link graph
  dbdet_edgel_link_graph& get_link_graph() { return edge_link_graph_; }

  //: return the curve fragment graph
  dbdet_curve_fragment_graph & get_curve_fragment_graph() { return curve_frag_graph_; }

  //********************************************************************//
  // User Friendly functions
  //********************************************************************//

  //: use the recommended sub-algorithms to extract final contour set (naive users should call this function)
  void extract_image_contours();

  //********************************************************************//
  // Functions to perform local grouping of edgemap_->edgels
  //********************************************************************//

  //: add a newly formed curvelet to all the relevant lists
  void add_a_curvelet(dbdet_curvelet* new_cvlet, unsigned ref_edge_ind, dbdet_curvelet* cvlet);

  //: form curvelets in a hierarchical fashion using combination rules
  void build_curvelets_using_combination_rules();
    void build_pairs();
    void build_triplets();
    void build_quadruplets();

    virtual void form_an_edgel_pair(dbdet_edgel* /*e1*/, dbdet_edgel* /*e2*/){}
    virtual void form_an_edgel_triplet(dbdet_curvelet* /*p1*/, dbdet_curvelet* /*p2*/){}
    virtual void form_an_edgel_quad(dbdet_curvelet* /*t1*/, dbdet_curvelet* /*t2*/){}

    //: add a newly formed edgel triplet to all relevant lists
    void add_a_triplet(dbdet_curvelet* trip, unsigned pos, dbdet_curvelet* p1, dbdet_curvelet* p2);
    //: add a newly formed edgel quadruplet to all relevant lists
    void add_a_quadruplet(dbdet_curvelet* quad, unsigned pos1, unsigned pos2, dbdet_curvelet* t1, dbdet_curvelet* t2);

    //: form curvelets in a hierarchical fashion (as high as the neighborhood allows)
    void build_curvelets_hierarchically();
    //: form an edgel grouping by adding a pair to the existing curvelet
    void form_an_edgel_grouping(dbdet_edgel* /*eA*/, dbdet_curvelet* /*cvlet*/, dbdet_curvelet* /*pair*/){}
    
    //: form curvelets around each edgel in a greedy fashion
    void build_curvelets_greedy(unsigned max_size_to_group, bool use_flag=false,  bool clear_existing=true, bool verbose=false);
    //: form curvelets around the given edgel in a greedy fashion
    virtual void build_curvelets_greedy_for_edge(dbdet_edgel* eA, unsigned max_size_to_group, 
        bool use_flag=false, bool forward=true,  bool centered=true, bool leading=true) = 0;

    //: form an edgel grouping from an ordered list of edgemap_->edgels
    virtual dbdet_curvelet* form_an_edgel_grouping(dbdet_edgel* ref_e, vcl_deque<dbdet_edgel*> &edgel_chain, 
        bool forward=true,  bool centered=true, bool leading=true) = 0;

    //: check to see if curvelets are balanced
    bool curvelet_is_balanced(dbdet_edgel* ref_e, vcl_deque<dbdet_edgel*> &edgel_chain);

    //: form the full curvelet map (curvelet map lists all the curvelets it participated in and not just the ones anchored to it)
    void form_full_cvlet_map();
  
    //: prune the curvelets that are below the quality threshold and hence considered spurious
    void prune_the_curvelets(double quality_threshold);
    void recompute_curvelet_quality();
    
    //: prune the curvelets with lengths (extent) smaller than the one specified
    void prune_curvelets_by_length(double length_threshold);
  
    //: prune the curvelets with gaps larger than the one specified
    void prune_curvelets_by_gaps(double gap_threshold);
  
    //: prne the curvelets that are not locally geometrically consistent (i.e., c1)
    void prune_curvelets_by_c1_condition();
  
    //: set the flags to detecmine if appearance is to be used for forming curvelets
    void set_appearance_usage(unsigned app_usage) { app_usage_ = app_usage; }
    void set_appearance_threshold(double app_thresh) { app_thresh_ = app_thresh; }
  
    //************************************************************************************//
    // Functions for forming a Link Graph from Curvelets perform local grouping of edgemap_->edgels
    //************************************************************************************//
  
    //: construct a simple link grpah by connectng edgels to all its neighbors
    void construct_naive_link_graph(double proximity_threshold=3.0, double affinity_threshold=10.0);
  
    //: form the link graph from the existing edgel groupings
    void construct_the_link_graph(unsigned min_group_size=4, int method=0);

    void use_all_curvelets() { use_anchored_curvelets_ = false;}
    void use_anchored_curvelets_only() { use_anchored_curvelets_ = true; }

    //: count the degree of overlap between two curvelets
    int count_degree_overlap(dbdet_curvelet* cvlet1, dbdet_curvelet* cvlet2, unsigned k1, unsigned k2);
    //: count the degree of overlap between pairs of curvelets at a link
    int count_degree_overlap(dbdet_link* link);

    //: form all appropriate links from a curvelet around an edgel
    void form_links_from_a_curvelet(dbdet_edgel* eA, dbdet_curvelet* cvlet, unsigned min_group_size=4, int method=0);

    //: check to see if link between two edges is reciprocal (i.e., x-A->b-y && p-a->B-q)
    bool link_is_reciprocal(dbdet_edgel* eA, dbdet_edgel* eB, unsigned min_group_size);
    bool link_is_reciprocal(dbdet_edgel* eA, dbdet_edgel* eB, dbdet_edgel* ref_e, unsigned min_group_size);

    //: check to see if link is reciprocal and supported by other edgemap_->edgels (i.e., x-A->b-y && x-a->B-y)
    bool link_is_supported(dbdet_edgel* eA, dbdet_edgel* eB, unsigned min_group_size);
    bool link_is_supported(dbdet_edgel* eX, dbdet_edgel* eA, dbdet_edgel* eB, dbdet_edgel* eY, dbdet_edgel* ref_e, 
                           unsigned min_group_size);

    //: check to see if this triplet is supported
    bool triplet_is_supported(dbdet_edgel* eX, dbdet_edgel* eA, dbdet_edgel* eY, 
                              unsigned min_group_size);

    //: prune the link graph by removing links from sub-standard cvlets
    void prune_the_link_graph();

    //: make the link graph bidirectionally consistent
    void make_link_graph_consistent();

    //: determine if this link is bidirectional
    bool link_bidirectional(dbdet_link* link);

  //********************************************************************//
  // Functions for extracting image contours
  //********************************************************************//

  //: clear all contours
  void clear_all_contours();

  //: extract image contours from the link graph by
  //  extracting regular contours in successive stages
  void extract_image_contours_from_the_link_graph(unsigned num_link_iters=1);

  //: extract the regular contours from the link graph
  void extract_regular_contours_from_the_link_graph();

    //: extract the one chains from a given link graph
    void extract_one_chains_from_the_link_graph(dbdet_edgel_link_graph& ELG);
  
    //: determine if this links from the link grpah is valid
    bool link_valid(dbdet_link* link);
    //: set the min_deg_to_link
    void set_min_deg_to_link(unsigned deg){ min_deg_to_link_ = deg; }
    //: get the min deg to link
    unsigned min_deg_to_link() { return min_deg_to_link_; }
  
  //********************************************************************//
  // Functions for fitting C^1 polyarcs to edgel chains
  //********************************************************************//

  //: fit C^1 polyarcs to all the unambiguous edgel chains in the image
  void fit_polyarcs_to_all_edgel_chains();

  //: attempt to fit polyarcs to edgel chains (knot-based algorithm)
  void fit_polyarc_to_chain(dbdet_edgel_chain* chain);


  //***********************************************************************//
  // Hypothesis trees and Hypothesis graphs
  //***********************************************************************//

  //Preprocessing of the link graph for the HTG algorithm

  //: seperate the link graph into two separate ones based on linking direction
  void separate_link_graphs_and_cvlet_maps();
    
    //: DFS search over the links to label the links and curvelets in a connected component scheme
    void DFS_CC_links(dbdet_link* link);
    //: BFS search over the links to label the links and curvelets in a connected component scheme
    void BFS_CC_links(dbdet_link* link);

  //: clear the HTG and related data structures
  void clear_HTG();

  //: initialize the hypothesis tree graph (basically mark the root nodes of each tree)
  void initialize_HTG();
    bool cvlet_has_bifurcations(dbdet_curvelet* cvlet, bool dir);

  //: propagate all the hypothesis trees
  void propagate_all_HTs();
  void propagate_HTs_N_steps(int N);

    void propagate_HT_from_the_next_leaf_node();
      vcl_set<dbdet_hyp_tree_node*> check_for_interaction(dbdet_hyp_tree_node* node);
      bool C1_CPL_interaction(dbdet_hyp_tree_node* node1, dbdet_hyp_tree_node* node2);

    //bool explore_continuations(vcl_queue<dbdet_hyp_tree_node*>& BFS_queue, dbdet_edgel* e);
    void explore_continuations(dbdet_hyp_tree_node* cur_node, dbdet_link* link);

      int dist_from_edge(dbdet_curvelet* cvlet, dbdet_edgel* e);
      bool cvlet_advances_curve(dbdet_curvelet* last_cvlet, dbdet_curvelet* cur_cvlet);
      bool cvlets_overlap_legally(dbdet_curvelet* last_cvlet, dbdet_curvelet* cur_cvlet);
      dbdet_curvelet* construct_chopped_cvlet(dbdet_curvelet* cvlet, dbdet_edgel* e);

      //: determine if C1 continuation is possible between a pair of curvelets
      dbdet_curve_model* C1_continuity_possible(dbdet_curvelet* cur_cvlet, dbdet_curvelet* next_cvlet);

      void label_edgels(dbdet_hyp_tree_node* node);
      void remove_labels(dbdet_hyp_tree_node* node);

      void remove_node_if_redundant(dbdet_hyp_tree_node* cur_node);

  //: perform gradient descent disambiguation of the HTG to segment the contours
  void disambiguate_the_HTG();

    void resolve_HTG_completion();
      void determine_optimal_HTG_completion_path(dbdet_HTG_link_path* link);
      void resolve_HTG_completion_path(dbdet_HTG_link_path* link);

    void resolve_HTG_competition();
      bool prune_HTs_of_claimed(dbdet_hyp_tree* HT);
        void delete_HT_subtree(dbdet_hyp_tree_node* cur_node);
      double determine_best_path(dbdet_hyp_tree* HT);
      void claim_best_path(dbdet_hyp_tree* HT);
        void claim_edgels(dbdet_hyp_tree_node* node);

    void resolve_HT(dbdet_hyp_tree* HT); 

    double compute_path_metric(vcl_vector<dbdet_curvelet*>& path);
    void back_propagate_solution(vcl_vector<dbdet_curvelet*>& path, vgl_point_2d<double> sol);
    
  void print_all_trees();

  //***********************************************************************//
  // Edgel based Hypothesis trees and CFTGs
  //***********************************************************************//

  //: construct a hyp tree from a given edgel 
  dbdet_EHT* construct_hyp_tree(dbdet_edgel* e);

  //: by Yuliang, a filter befor building hypothesis tree, rule out some strang pathe from regular contours, and merge some contours, not deal with junction
  void regular_contour_filter();

  //: construct all possible EHTs from the terminal nodes and find legal contour paths
  void construct_all_path_from_EHTs();

  //: Construction of Hypothesis Tree ::// Naman Kumar
  void Construct_Hypothesis_Tree();

  //:Disambiguating the Hypothesis Tree  ::// Naman Kumar
  void Disambiguation();

  //: Minor post Process to prune some extra small part of contours ::// Naman Kumar
  void Post_Process();

  //: perform a geometric consistency check to determine whether a given temp path is valid
  bool is_EHT_path_legal(vcl_vector<dbdet_edgel*>& edgel_chain);

  //: compute a simple path metric based on the chain and its neighboring support chains
  double compute_path_metric2(vcl_vector<dbdet_edgel*>& Pchain, 
                             vcl_vector<dbdet_edgel*>& Tchain, 
                             vcl_vector<dbdet_edgel*>& Cchain);

  //: disambiguate the CFG, basically to produce a disjoint set
  void disambiguate_the_CFTG();

  void merge_extreme_short_curve_frags();

  //: correct the CFG topology to produce a disjoint set
  void correct_CFG_topology();

  //: Construct Hypothesis Tree to improve Contour Extraction
  void Hypothesis_tree_construction();

  //***********************************************************************//
  // Hybrid method: Perform geometric consistency check on curvefragments
  //***********************************************************************//

  //: construct a mapping between edgel ids and curve ids
  void compile_edge_to_contour_mapping();

  void form_curvelets_from_contours(bool clear_existing=true);

  //: attempt to form curvelets from the traced contour fragments
  void form_curvelets_from_contours(unsigned max_size_to_group);

  //: Break contours that do not agree with the curvelets
  void post_process_to_break_contours();

  //: connect pieces of contours that are connected but separated
  void post_process_to_link_contour_fragments();

  //********************************************************************//
  // Miscellaneous functions
  //********************************************************************//

  //: method to look at the spread of differential estimates
  void determine_accuracy_of_measurements();

  //: for debug
  void report_stats();

  //: evauate the qualities of curvelets using various functions
  void evaluate_curvelet_quality(int method);

protected:

  dbdet_edgemap_sptr edgemap_; ///< the edgemap to link
  dbdet_curvelet_map& curvelet_map_;             ///< The curvelet map (CM)
  dbdet_edgel_link_graph& edge_link_graph_;      ///< The edge link graph (ELG)
  dbdet_curve_fragment_graph& curve_frag_graph_; ///< The curve fragment graph (CFG)

  unsigned nrows_, ncols_;     ///< size of the edgemap

  //curvelet formation parameters
  unsigned app_usage_;
  double app_thresh_; ///< appearance threshold
  double rad_;    ///< radius of the grouping neighborhood around each edgel
  double dtheta_; ///< assumed uncertainty in orientation
  double dpos_;   ///< assumed uncertainty in position
  bool badap_uncer_; ///< use the uncertainty stored in the edgels
  double token_len_; ///< Length of the edgel token (puts bounds on legal curvature)
  double max_k_;  ///< maximum curvature of curves in the curve bundle
  double max_gamma_; ///< maximum curvature derivative of curves in the curve bundle

  unsigned nrad_; ///< the range of edgel cells to search (this should contain all edgemap_->edgels within rad_) 
  double gap_; ///<Distance between two consecutive edges in a link
  unsigned maxN_; ///< largest curvelet size to form
  bool centered_; ///< curvelets centered on the anchor edgel
  bool bidir_;    ///< curvelets in both direction

  //linking parameters
  bool use_anchored_curvelets_; ///< the curvelet set to use for linking
  unsigned min_deg_to_link_; ///< minimum degree of a link before it is linked 

  //for the connected components algo to separate the link graph and the cvlet map
  vcl_set<dbdet_curvelet*> cv_set1;
  vcl_map<vcl_pair<int, int>, dbdet_link*> link_map;
  vcl_queue<dbdet_link*> BFS_links_queue;
  vcl_map<dbdet_link*, dbdet_link*> link_pairs;

public: //temp

  //for the hybrid algorithm
  bool use_hybrid_; 
  vcl_vector<int> cId_; ///< id of the contour that each edgel belongs to

  //for the DHT algorithm
  dbdet_HT_graph HTG;             ///< The global HTG
  vcl_queue<dbdet_hyp_tree_node*> BFS_queue_global; ///< BFS queue for propagating the HTs simultaneously

  dbdet_edgel_labels ELs;         ///< edgel labels for contour ownership

  bool DHT_mode_; ///< there are some steps that are experimental and only applies to DHT linking

  bool propagate_constraints;
};

//: data structure to hold pairwise hypothesis around an edgel
struct sel_hyp {
  bool ref_first;           ///< the direction of groouping    
  double d;                 ///< distance between them
  bool flag;                ///< if this hypothesis has been considered before
  dbdet_edgel* eN;          ///< neighboring edgel
  dbdet_curve_model* cm;    ///< the curve bundle of the appropriate model
};

inline bool comp_dist_hyps_less(const sel_hyp &h1, const sel_hyp &h2)
{
  if (h1.d < h2.d)
    return true;
  else
    return false;
}

#endif // dbdet_sel_base_h
