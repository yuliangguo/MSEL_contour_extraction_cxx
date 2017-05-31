// This is brcv/seg/dbdet/pro/dbdet_sel_process.cxx
// modified by Yuliang Jul/29/2010  can choose not to form link graph
//:
// \file

#include "dbdet_sel_process.h"

#include "dbdet_edgemap_storage.h"
#include "dbdet_edgemap_storage_sptr.h"
#include "dbdet_edgemap_sptr.h"

#include "dbdet_sel_storage.h"
#include "dbdet_sel_storage_sptr.h"
#include "dbdet_sel_sptr.h"
#include "dbdet_sel.h"
#include "dbdet_curve_model.h"

#include <vcl_vector.h>
#include <vcl_string.h>
#include <vul/vul_timer.h>
#include <vil/vil_image_resource.h>

//: Constructor
dbdet_sel_process::dbdet_sel_process()
{
  vcl_vector<vcl_string> curve_model_choices;
  curve_model_choices.push_back("Simple Linear model");                   //0
  curve_model_choices.push_back("Linear model");                          //1
  curve_model_choices.push_back("Circular Arc w/o perturbations");        //2
  curve_model_choices.push_back("Circular Arc w k classes");              //3
  curve_model_choices.push_back("Circular Arc w  perturbations");         //4
  curve_model_choices.push_back("Circular Arc 3d Bundle");                //5
  curve_model_choices.push_back("ES w/o perturbations");        //6
  curve_model_choices.push_back("ES w   perturbations");        //7

  //grouping_algo_choices
  vcl_vector<vcl_string> grouping_algo_choices;
  grouping_algo_choices.push_back("Combinatorial grouping");               //0
  grouping_algo_choices.push_back("Hierarchical grouping");                //1
  grouping_algo_choices.push_back("Greedy Local grouping");                //2
  grouping_algo_choices.push_back("Very Greedy Local grouping");           //3

  //Curvelet type choices
  vcl_vector<vcl_string> curvelet_type_choices;
  curvelet_type_choices.push_back("Anchor Centered");                      //0
  curvelet_type_choices.push_back("Anchor Centered/Bidirectional");        //1
  curvelet_type_choices.push_back("Anchor Leading/Bidirectional");         //2
  curvelet_type_choices.push_back("ENO Style around Anchor");              //3

  //Appearance_choices
  vcl_vector<vcl_string> appearance_usage_choices;
  appearance_usage_choices.push_back("Do Not Use Appearance");             //0
  appearance_usage_choices.push_back("Use Local Comparison");              //1
  appearance_usage_choices.push_back("Compare against Ref.");              //2

  //LinkGraph_choices
  vcl_vector<vcl_string> linkgraph_choices;
  linkgraph_choices.push_back("all links");                                //0
  linkgraph_choices.push_back("immediate links only");                     //1
  linkgraph_choices.push_back("immediate reciprocal links");               //2
  linkgraph_choices.push_back("immediate reciprocal links with support");  //3
  linkgraph_choices.push_back("triplets with support");                    //4

  //Linking algo choices
  vcl_vector<vcl_string> linking_algo_choices;
  linking_algo_choices.push_back("Do not Link");                           //0
  linking_algo_choices.push_back("From the Link Graph");                   //1
  linking_algo_choices.push_back("Regular Contours");                      //2

  if ( 
      //grouping parameters
      !parameters()->add( "Radius of Neighborhood" , "-nrad" , 3.5 ) ||
      //By Naman Kumar
      !parameters()->add( "Maximum Pixel Distance to complete" , "-gap" , 3.0 ) ||
      !parameters()->add( "Get Uncertainty from edges" , "-badap_uncer" , true ) ||
      !parameters()->add( "  - Position uncertainty" , "-dx" , 0.4 ) ||
      !parameters()->add( "  - Orientation uncertainty(Deg)" , "-dt" , 15.0 ) ||
      
      //curve model
      !parameters()->add( "Curve Model"   , "-curve_model" , curve_model_choices, 5) ||

      //curve model parameters
      !parameters()->add( "  - Token Length" , "-token_len" , 1.0 ) ||
      //!parameters()->add( "  - Maximum Curvature" , "-max_k" , 0.2 ) ||
      //!parameters()->add( "  - Maximum Curvature Derivative" , "-max_gamma" , 0.05 ) ||

      //grouping algorithm
      !parameters()->add( "Local Grouping Algorithm"   , "-grouping_algo" , grouping_algo_choices, 2) ||
      !parameters()->add( "  - Curvelet Type"   , "-cvlet_type" , curvelet_type_choices, 0) ||
      !parameters()->add( "  - Appearance Usage"   , "-app_usage" , appearance_usage_choices, 0) ||
      !parameters()->add( "  - App Threshold" , "-app_thresh" , 0.2 ) ||

      !parameters()->add( "  - Maximum # of edgels to group" , "-max_size_to_group", (unsigned) 7 ) ||

      !parameters()->add( "Form Complete Curvelet Map" , "-bFormCompleteCvletMap", false ) ||

      !parameters()->add( "Form Link Graph", "-bFormLinkGraph", true ) ||

      //link graph formation
      //!parameters()->add( "Use All Curvelets" , "-b_use_all_cvlets", false ) ||

      //!parameters()->add( "Form Link Graph Using"   , "-linkgraph_algo" , linkgraph_choices, 0) ||
      //!parameters()->add( "Minimum curvelet size to link" , "-min_size_to_link", (unsigned) 4 ) ||

      //link graph --> image contours
      !parameters()->add( "Extract Image Contours "   , "-linking_algo" , linking_algo_choices, 0 ) ||
      
      //linking parameter for contour extraction from the link graph
      !parameters()->add( "  - Number of Linking iterations" , "-num_link_iters", (unsigned) 7 ) ||
      
      //Get Final contours in one step
      !parameters()->add( "Proceed with Hypothesis Graph", "-bGetfinalcontours", true ) ||

	  //merge frag candidates after hypothesis
      !parameters()->add( "Merge curve fragments candidates", "-bmergefrags", true )
	)
  {
    vcl_cerr << "ERROR: Adding parameters in " __FILE__ << vcl_endl;
  }
}


//: Destructor
dbdet_sel_process::~dbdet_sel_process()
{
}


//: Clone the process
bpro1_process*
dbdet_sel_process::clone() const
{
  return new dbdet_sel_process(*this);
}


//: Return the name of this process
vcl_string
dbdet_sel_process::name()
{
  return "Symbolic Edge Linker";
}


//: Return the number of input frame for this process
int
dbdet_sel_process::input_frames()
{
  return 1;
}


//: Return the number of output frames for this process
int
dbdet_sel_process::output_frames()
{
  return 1;
}


//: Provide a vector of required input types
vcl_vector< vcl_string > dbdet_sel_process::get_input_type()
{
  vcl_vector< vcl_string > to_return;
  to_return.push_back( "edge_map" );

  return to_return;
}


//: Provide a vector of output types
vcl_vector< vcl_string > dbdet_sel_process::get_output_type()
{
  vcl_vector<vcl_string > to_return;
  //output the sel storage class
  to_return.push_back( "sel" );
  return to_return;
}


//: Execute the process
bool dbdet_sel_process::execute()
{
  if ( input_data_.size() != 1 ){
    vcl_cout << "In dbdet_sel_process::execute() - not exactly one input \n";
    return false;
  }
  clear_output();

  //get the parameters
  get_parameters();

  //get the input storage class
  dbdet_edgemap_storage_sptr input_edgemap;
  input_edgemap.vertical_cast(input_data_[0][0]);
  dbdet_edgemap_sptr EM = input_edgemap->get_edgemap();

  // create the sel storage class
  dbdet_sel_storage_sptr output_sel = dbdet_sel_storage_new();
  output_sel->set_EM(EM);
  //get pointers to the data structures in it
  dbdet_curvelet_map& CM = output_sel->CM();
  dbdet_edgel_link_graph& ELG = output_sel->ELG();
  dbdet_curve_fragment_graph &CFG = output_sel->CFG();

  //different types of linkers depending on the curve model
  typedef dbdet_sel<dbdet_simple_linear_curve_model> dbdet_sel_simple_linear;
  typedef dbdet_sel<dbdet_linear_curve_model> dbdet_sel_linear;
  typedef dbdet_sel<dbdet_CC_curve_model> dbdet_sel_CC;
  typedef dbdet_sel<dbdet_CC_curve_model_new> dbdet_sel_CC_new;
  typedef dbdet_sel<dbdet_CC_curve_model_perturbed> dbdet_sel_CC_perturbed;
  typedef dbdet_sel<dbdet_CC_curve_model_3d> dbdet_sel_CC_3d;
  typedef dbdet_sel<dbdet_ES_curve_model> dbdet_sel_ES;
  typedef dbdet_sel<dbdet_ES_curve_model_perturbed> dbdet_sel_ES_perturbed;

  //start the timer
  vul_timer t;

  //construct the linker
  dbdet_sel_sptr edge_linker;

  //The curvelet formation parameters
  dbdet_curvelet_params cvlet_params(dbdet_curve_model::CC, 
                                     nrad, gap, dt, dx, badap_uncer, 
                                     token_len, max_k, max_gamma,
                                     bCentered_grouping,
                                     bBidirectional_grouping);

  switch (curve_model_type)
  {
  case 0: //simple linear_model
    cvlet_params.C_type = dbdet_curve_model::LINEAR;
    edge_linker = new dbdet_sel_simple_linear(EM, CM, ELG, CFG, cvlet_params);
    break;
  case 1: //linear_model
    cvlet_params.C_type = dbdet_curve_model::LINEAR;
    edge_linker = new dbdet_sel_linear(EM, CM, ELG, CFG, cvlet_params);
    break;
  case 2: //CC_model
    cvlet_params.C_type = dbdet_curve_model::CC;
    edge_linker = new dbdet_sel_CC(EM, CM, ELG, CFG, cvlet_params);
    break;
  case 3: //CC_model new
    cvlet_params.C_type = dbdet_curve_model::CC2;
    edge_linker = new dbdet_sel_CC_new(EM, CM, ELG, CFG, cvlet_params);
    break;
  case 4: //CC_model with discrete perturbations
    cvlet_params.C_type = dbdet_curve_model::CC;
    edge_linker = new dbdet_sel_CC_perturbed(EM, CM, ELG, CFG, cvlet_params);
    break;
  case 5: //CC_model 3d bundle
    cvlet_params.C_type = dbdet_curve_model::CC3d;
    edge_linker = new dbdet_sel_CC_3d(EM, CM, ELG, CFG, cvlet_params);
    break;
  case 6: //ES_model
    cvlet_params.C_type = dbdet_curve_model::ES;
    edge_linker = new dbdet_sel_ES(EM, CM, ELG, CFG, cvlet_params);
    break;
  case 7: //ES_model with discrete perturbations
    cvlet_params.C_type = dbdet_curve_model::ES;
    edge_linker = new dbdet_sel_ES_perturbed(EM, CM, ELG, CFG, cvlet_params);
    break;
  }
  
  //set appearance usage flags
  edge_linker->set_appearance_usage(app_usage);
  edge_linker->set_appearance_threshold(app_thresh);

  //perform local edgel grouping
  switch (grouping_algo)
  {
  case 0: //combinatorial grouping
    edge_linker->build_curvelets_using_combination_rules();
    break;
  case 1: //hierarchical grouping (breadth-first grouping)
    edge_linker->build_curvelets_hierarchically();
    break;
  case 2: //greedy (depth first grouping)
    edge_linker->build_curvelets_greedy(max_size_to_group, false);
    break;
  case 3: //extra greedy (depth first grouping)
    edge_linker->build_curvelets_greedy(max_size_to_group, true);
    break;
  }
  
  if (bFormCompleteCvletMap)
    edge_linker->form_full_cvlet_map();
  
  double group_time = t.real() / 1000.0;
  t.mark();
  vcl_cout << "Time taken to form groups and cunstruct curvelet map: " << group_time << " sec" << vcl_endl;

  if (bFormLinkGraph){
  //form a link graph
  if (b_use_all_cvlets)
    edge_linker->use_all_curvelets();
  else
    edge_linker->use_anchored_curvelets_only();

  //form the link graph
  edge_linker->construct_the_link_graph(min_size_to_link, linkgraph_algo);
  
  //extract contours
  switch (linking_algo) {
    case 1: // iteratively
      edge_linker->extract_image_contours_from_the_link_graph(num_link_iters);
      break;
    case 2: // regular contours
      edge_linker->extract_regular_contours_from_the_link_graph();
      break;
  }
// By Yuliang Guo, Oct, 2010
      edge_linker->extract_regular_contours_from_the_link_graph();
  if(bGetfinalcontours){  
      edge_linker->Construct_Hypothesis_Tree();
      edge_linker->Disambiguation();
      if(bmergefrags){
		  edge_linker->correct_CFG_topology(); 
		  //edge_linker->Post_Process();
	  }
	  else
	  {
		 edge_linker->merge_extreme_short_curve_frags();
		 //edge_linker->Post_Process();
	  }
      //edge_linker->correct_CFG_topology();
  }
  double link_time = t.real() / 1000.0;
  vcl_cout << "Time taken to link: " << link_time << " sec" << vcl_endl;
}
  //report stats
  //edge_linker->report_stats();
  //edge_linker->determine_accuracy_of_measurements();

  // output the sel storage class
  output_data_[0].push_back(output_sel);

  return true;
}

bool
dbdet_sel_process::finish()
{
  return true;
}

void
dbdet_sel_process::get_parameters()
{
  //default parameter values
  token_len = 1.0;
  max_k = 0.2;
  max_gamma = 0.05;
  bCentered_grouping = true;

  linkgraph_algo = 0;
  min_size_to_link = 4;
  b_use_all_cvlets = false;
  linking_algo = 0;
  app_usage = 0;
  app_thresh = 0.2;

  //grouping parameters
  parameters()->get_value( "-nrad", nrad);
  parameters()->get_value( "-gap", gap);
  parameters()->get_value( "-badap_uncer", badap_uncer);
  parameters()->get_value( "-dx", dx);
  parameters()->get_value( "-dt", dt);

  //curve model
  parameters()->get_value( "-curve_model" , curve_model_type);

  //curve model parameters
  parameters()->get_value( "-token_len" , token_len );
  //parameters()->get_value( "-max_k" , max_k );
  //parameters()->get_value( "-max_gamma" , max_gamma );

  //grouping algorithm 
  parameters()->get_value( "-grouping_algo" , grouping_algo);
  parameters()->get_value( "-cvlet_type" , cvlet_type);

  switch(cvlet_type) //set the grouping flags from the choice of cvlet type
  {
  case 0: //Anchor Centered
    bCentered_grouping = true;
    bBidirectional_grouping = false;
    break;
  case 1: //Anchor Centered/Bidirectional
    bCentered_grouping = true;
    bBidirectional_grouping = true;
    break;
  case 2: //Anchor Leading/Bidirectional
    bCentered_grouping = false;
    bBidirectional_grouping = true;
    break;
  case 3: //ENO Style around Anchor
    bCentered_grouping = false;
    bBidirectional_grouping = false;
    break;
  }
  parameters()->get_value( "-max_size_to_group", max_size_to_group);

  //appearance usage options
  parameters()->get_value( "-app_usage" , app_usage);
  parameters()->get_value( "-app_thresh" , app_thresh );
  
  parameters()->get_value( "-bFormCompleteCvletMap" , bFormCompleteCvletMap );

  parameters()->get_value( "-bFormLinkGraph" , bFormLinkGraph );
  //algorithm to form the link graph
  //parameters()->get_value( "-linkgraph_algo" , linkgraph_algo);
  //parameters()->get_value( "-min_size_to_link", min_size_to_link);
  //parameters()->get_value( "-b_use_all_cvlets", b_use_all_cvlets);

  //extract image contours from the link graph
  parameters()->get_value( "-linking_algo" , linking_algo);
  parameters()->get_value( "-num_link_iters" , num_link_iters);
  parameters()->get_value( "-bGetfinalcontours" , bGetfinalcontours );
  parameters()->get_value( "-bmergefrags" , bmergefrags );
  
}

