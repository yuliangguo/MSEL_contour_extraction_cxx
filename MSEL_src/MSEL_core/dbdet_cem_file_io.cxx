#include "dbdet_cem_file_io.h"

#include "dbdet_edgemap.h"
#include "dbdet_curve_fragment_graph.h"
#include "dbdet_postprocess_contours.h"
#include <vcl_cstring.h>
#include <vcl_cassert.h>

// -----------------------------------------------------------------------------
//: Save the contour fragment graph as a .cem file

static bool sort_cmp(dbdet_edgel * e1, dbdet_edgel * e2)
{
  return e1->id < e2->id;
}

bool dbdet_save_cem(vcl_string filename, dbdet_edgemap_sptr EM, dbdet_curve_fragment_graph& CFG)
{
  //1)If file open fails, return.
  vcl_ofstream outfp(filename.c_str(), vcl_ios::out);

  if (!outfp){
    vcl_cout << " Error opening file  " << filename.c_str() << vcl_endl;
    return false;
  }

/*
  // Compile edge information first
  vcl_vector<int> emap(EM->num_edgels(), -2); //-2 indicates unused

  dbdet_edgel_chain_list_iter cit = CFG.frags.begin();
  for (; cit!=CFG.frags.end(); cit++){
    dbdet_edgel_chain* chain = (*cit);

    for (unsigned i=0; i<chain->edgels.size(); i++)  
      emap[chain->edgels[i]->id] = -1; //-1 indicates used
  }

  //compile new edge list (only keep the used ones in the new list)
  vcl_vector<int> new_emap;
  for (unsigned i=0; i<EM->edgels.size(); i++){
    if (emap[i]==-1){
      emap[i] = new_emap.size(); //save the forward mapping
      new_emap.push_back(i);     //save the reverse mapping
    }
  }
*/
  // change by Yuliang, Oct, 2014, save the full edge map istead of the used edges
  //vcl_vector<int> emap(EM->num_edgels(), 0); //-2 indicates unused

//compile new edge list (only keep the used ones in the new list)
 /* vcl_vector<int> new_emap;
  for (unsigned i=0; i<EM->edgels.size(); i++){
//    if (emap[i]==-1){
      emap[i] = new_emap.size(); //save the forward mapping
      new_emap.push_back(i);     //save the reverse mapping
    }
*/

  // output header information
  outfp << ".CEM v2.0 " << vcl_endl;
  outfp << "size=[" << EM->width() << " " << EM->height() << "]" << vcl_endl;

  //output the edgemap section
  outfp << "[Edgemap]" << vcl_endl;
  outfp << "count=" << EM->edgels.size() << vcl_endl;
  //outfp << "# (x, y) dir conf d2f" << vcl_endl;

  for (unsigned i=0; i<EM->edgels.size(); i++){
    dbdet_edgel* e = EM->edgels[i];
    assert(i == e->id);
    outfp << "(" << e->pt.x() << ", " << e->pt.y() << ")\t" << e->tangent << "\t" << e->strength << "\t" << e->deriv << vcl_endl; 
  }

  //output the contours section
  outfp << "[Contours]" << vcl_endl;
  outfp << "count=" << CFG.frags.size() << vcl_endl;

  dbdet_edgel_chain_list_iter cit = CFG.frags.begin();
  for (; cit!=CFG.frags.end(); cit++){
    dbdet_edgel_chain* chain = (*cit);

    outfp << "[";
    for (unsigned i=0; i<chain->edgels.size(); i++)
      outfp << chain->edgels[i]->id << " ";
    outfp << "]" << vcl_endl;
  }

  //next output some properties of the contours that could be useful for evaluation
  outfp << "[Contour Properties]" << vcl_endl;
  outfp << "# <len> <avg. str> <mean con> <Lstd> <Rstd> <avg. d2f> <avg. k> <max k>" << vcl_endl;

  cit = CFG.frags.begin();
  for (; cit!=CFG.frags.end(); cit++){
    dbdet_edgel_chain* chain = (*cit);

    //compute various properties
    outfp << chain->edgels.size() << " ";
    outfp << avg_strength(chain) << " ";
    outfp << mean_contrast(chain) << " ";
    outfp << left_app_std(chain) << " ";
    outfp << right_app_std(chain) << " ";
    outfp << avg_d2f(chain) << " ";
    outfp << avg_curvature(chain) << " ";
    outfp << max_curvature(chain) << " ";

    outfp << vcl_endl;
  }


  //close file
  outfp.close();
  
  return true;
}


// -----------------------------------------------------------------------------
//: load the older version of the cem file (for backward compatibility)
dbdet_edgemap_sptr dbdet_load_cem_v1(vcl_ifstream &infp, dbdet_curve_fragment_graph& CFG, int width, int height, bool convert_degrees_to_radians)
{
  char lineBuffer[1024];
  int numContours, numTotalEdges, numEdges;
  int ix, iy;
  double x, y;
  double idir, iconf, dir, conf;

  //construct the edgemap
  dbdet_edgemap_sptr edgemap = new dbdet_edgemap(width, height); //CEM file does not have size info
  vcl_vector<dbdet_edgel_chain*> chains;

  //2)Read in each line
  while (infp.getline(lineBuffer,1024)) {

    //ignore comment lines and empty lines
    if (vcl_strlen(lineBuffer)<2 || lineBuffer[0]=='#')
      continue;

    //read the line with the contour count info
    if (!vcl_strncmp(lineBuffer, "CONTOUR_COUNT=", sizeof("CONTOUR_COUNT=")-1)){
      sscanf(lineBuffer,"CONTOUR_COUNT=%d",&(numContours));
      //vcl_cout << numContours << vcl_endl;
      continue;
    }

    //read the line with the edge count info
    if (!vcl_strncmp(lineBuffer, "TOTAL_EDGE_COUNT=", sizeof("TOTAL_EDGE_COUNT=")-1)){
      sscanf(lineBuffer,"TOTAL_EDGE_COUNT=%d",&(numTotalEdges));
      //vcl_cout << numTotalEdges << vcl_endl;
      continue;
    }

    //read the beginning of a contour block
    if (!vcl_strncmp(lineBuffer, "[BEGIN CONTOUR]", sizeof("[BEGIN CONTOUR]")-1))
    {
      infp.getline(lineBuffer,1024);

      sscanf(lineBuffer,"EDGE_COUNT=%d",&(numEdges));
      //vcl_cout << numEdges << vcl_endl;

      dbdet_edgel_chain* chain = new dbdet_edgel_chain();

      for (int j=0; j< numEdges; j++){
        //the rest should have data that goes into the current contour
        infp.getline(lineBuffer,1024);

        sscanf(lineBuffer," [%d, %d]\t%lf\t%lf\t[%lf, %lf]\t%lf\t%lf",&(ix), &(iy),
              &(idir), &(iconf), &(x), &(y), &(dir), &(conf));

        if(convert_degrees_to_radians)
        {            //VJ's current CEM is in degrees rather than radians so need to convert
            dir += 90;
            dir *= vnl_math::pi/180;
        }

        dbdet_edgel* e = new dbdet_edgel(vgl_point_2d<double>(x,y), dir, conf);
        e->id = edgemap->edgels.size();
        edgemap->insert(e);

        chain->push_back(e);
      }

      //go to the end of the block
      infp.getline(lineBuffer,1024);
      while (vcl_strncmp(lineBuffer, "[END CONTOUR]", sizeof(" [END CONTOUR]")-1))
        infp.getline(lineBuffer,1024);

      //save it in the chains list for now
      chains.push_back(chain);
 
    }
  }
  infp.close();

  CFG.clear();
  CFG.resize(edgemap->num_edgels());
  for (unsigned i=0; i<chains.size(); i++)
    CFG.insert_fragment(chains[i]);

  return edgemap;
}




// -----------------------------------------------------------------------------
//: Loads an ascii file containing a graph of edgel chains (the contour fragment graph)
dbdet_edgemap_sptr dbdet_load_cem(vcl_string filename, dbdet_curve_fragment_graph& CFG)
{
  ////
  char lineBuffer[1024];
  int version=1;
  dbdet_edgemap_sptr edgemap;

  //1)If file open fails, return.
  vcl_ifstream infp(filename.c_str(), vcl_ios::in);

  if (!infp){
    vcl_cout << " Error opening file  " << filename.c_str() << vcl_endl;
    return NULL;
  }

  char * cur_locale, * dup_locale;
  cur_locale = setlocale(LC_NUMERIC, NULL);
  dup_locale = strdup(cur_locale);
  setlocale(LC_NUMERIC, "C");

  //determine the version of this file
  infp.getline(lineBuffer,1024); //read in the first line
  if (!vcl_strncmp(lineBuffer, ".CEM v2.0", sizeof(".CEM v2.0")-1))
  {
    version =2;
    edgemap = dbdet_load_cem_v2(infp, CFG);
    vcl_cout << "Loaded: " << filename.c_str() << ".\n";
  }
  else {
    edgemap = dbdet_load_cem_v1(infp, CFG);
    vcl_cout << "Loaded: " << filename.c_str() << ".\n";  
  }

  setlocale(LC_NUMERIC, dup_locale);
  free(dup_locale);

  return edgemap; 
}







// -----------------------------------------------------------------------------
//: load cem file version 2
dbdet_edgemap_sptr dbdet_load_cem_v2(vcl_ifstream &infp, dbdet_curve_fragment_graph& CFG)
{
  //
  char lineBuffer[1024];
  int width, height;
  double x, y, dir, conf, d2f;
  int num_contours = 0;

  dbdet_edgemap_sptr edgemap = 0;

  // Read in each line
  while (infp.getline(lineBuffer,1024)) {

    //ignore comment lines and empty lines
    if (vcl_strlen(lineBuffer)<2 || lineBuffer[0]=='#')
      continue;

    //read the line with the edgemap size info
    if (!vcl_strncmp(lineBuffer, "size=", sizeof("size=")-1)){
      sscanf(lineBuffer,"size=[%d %d]",&width, &height);
      
      //construct the edgemap
      edgemap = new dbdet_edgemap(width, height);

      continue;
    }

    //read the edgemap block
    if (!vcl_strncmp(lineBuffer, "[Edgemap]", sizeof("[Edgemap]")-1))
    {
      //read the next line with the edge count
      infp.getline(lineBuffer,1024);

      int num_edges = 0;
      sscanf(lineBuffer,"count=%d", &num_edges);

      //read in all the edges
      for (int i=0; i<num_edges; i++){
        infp.getline(lineBuffer,1024);

        if (sscanf(lineBuffer,"(%100lf, %100lf)\t%100lf\t%100lf\t%100lf",&x, &y, &dir, &conf, &d2f) != 5)
          std::cerr << "dbdet_load_cem_v2: input error for edgel #" << i << std::endl;

        //construct a new edgel and add it to the edgemap
        dbdet_edgel* e = new dbdet_edgel(vgl_point_2d<double>(x,y), dir, conf, d2f);
        e->id = edgemap->edgels.size();
        edgemap->insert(e);
      }
      continue;
    }

    //read the contours block
    if (!vcl_strncmp(lineBuffer, "[Contours]", sizeof("[Contours]")-1))
    {
      //read the next line with the edge count
      infp.getline(lineBuffer,1024);

      sscanf(lineBuffer,"count=%d", &num_contours);

      CFG.clear();
      CFG.resize(edgemap->num_edgels());

      //read in all the contours
      for (int i=0; i<num_contours; i++)
      {
        dbdet_edgel_chain* chain = new dbdet_edgel_chain();
        
        //we don't know how many edgels there are in this contour,
        //so we need to read in one by one
        unsigned char dummy;
        int e_id;
        infp >> dummy; //[
        infp >> e_id;
        
        while (!infp.fail()){ //read in all the ids until the end bracket is reached
          chain->push_back(edgemap->edgels[e_id]);
          infp >> e_id;
        }
        infp.clear(); //clear the fail bit

        //add this contour fragment to the graph
        CFG.insert_fragment(chain);

        //read in the line end character so that we can start at the 
        //beginning of the new line on the next iteration
        infp.getline(lineBuffer, 20);
      }
    }

    //read the contour properties block
    if (!vcl_strncmp(lineBuffer, "[Contour Properties]", sizeof("[Contour Properties]")-1))
    {
      //read the next line with the comment
      infp.getline(lineBuffer,1024);

      //read all the lines with the properties (but ignore them)
      for (int i=0; i<num_contours; i++)
        infp.getline(lineBuffer,1024);
    }
  }
  infp.close();
  return edgemap;
}





