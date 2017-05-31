// This is brcv/seg/dbdet/edge/dbdet_edgemap.h
#ifndef dbdet_edgemap_h
#define dbdet_edgemap_h
//:
//\file
//\brief The edge map class
//\author Amir Tamrakar
//\date 09/09/06
//
//\verbatim
//  Modifications
//\endverbatim

#include <vcl_vector.h>
#include <vcl_cmath.h>

#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_array_2d.h>
#include <vgl/vgl_point_2d.h>
#include "dbdet_edgel.h"
#include "dbdet_sel_utils.h"
#include <vnl/vnl_math.h>
#include <vcl_cmath.h>

//: very light class to store edgels
class dbdet_edge
{
public:
  vgl_point_2d<double> pt; ///< the location of the edgel
  double tangent;          ///< the orientation of the edgel
  double strength;         ///< the strength of the edgel (typically gradient magnitude)
  int id;                  ///< unique id

  dbdet_edge(vgl_point_2d<double> new_pt, double tan, double edge_strength=10.0):
    pt(new_pt), tangent(tan), strength(edge_strength) {}
  ~dbdet_edge(){}
};

typedef vbl_array_2d<vcl_vector<dbdet_edgel*> >::iterator dbdet_edgemap_iter;
typedef vbl_array_2d<vcl_vector<dbdet_edgel*> >::const_iterator dbdet_edgemap_const_iter;

//: A bucketing structure to hold the edgel tokens. Currently this is implemented
//  simply as a 2d array of tokens the same size as the image that gave rise to it.
//  For expediting linking from this data structure, I am using the dbdet_edgel class
//  to store the edgel tokens instead of the dbdet_edge class which is considerably lighter
//
class dbdet_edgemap : public vbl_ref_count
{
public:
  
  //: retinotopic map of edgels
  vbl_array_2d<vcl_vector<dbdet_edgel*> > edge_cells; 
  
  //: local list of edgels for easier traversal
  vcl_vector<dbdet_edgel*> edgels; 

  //: edgel occupancy map (redundant structure)
  vbl_array_2d<bool> occupancy;

  //: constructor
  dbdet_edgemap(int width, int height) : edgels(0) { edge_cells.resize(height, width); }

  //: constructor2
  dbdet_edgemap(int width, int height, vcl_vector<dbdet_edgel*>& edgels) : edgels(0)
  { 
    edge_cells.resize(height, width); 
    for (unsigned i=0; i<edgels.size(); i++)
      insert(edgels[i]);
  }

  //: destructor
  ~dbdet_edgemap()
  {
    //go over each cell and delete the edgels
    dbdet_edgemap_const_iter it = edge_cells.begin();
    for (; it!=edge_cells.end(); it++)
      for (unsigned j=0; j<(*it).size(); j++)
        delete (*it)[j];

    edge_cells.clear();

    //also clear the list of edgels
    edgels.clear();
  }

  //Access functions
  unsigned width() const { return edge_cells.cols(); }
  unsigned height() const { return edge_cells.rows(); }
  unsigned ncols() const  { return edge_cells.cols(); }
  unsigned nrows() const  { return edge_cells.rows(); }
  unsigned num_edgels() const { return edgels.size(); } ///< number of edgels in the edgemap

  //: read only access
  const vcl_vector<dbdet_edgel*>& cell(int x, int y) const { return edge_cells(y, x); }

  //: put an edgel into the edgemap at the prescribed cell
  void insert(dbdet_edgel* e, int xx, int yy)
  {
    edge_cells(yy, xx).push_back(e);
    e->id = edgels.size(); //assign unique id
    e->gpt = vgl_point_2d<int>(xx, yy); //record grid location

    edgels.push_back(e);
  }

  //: put an edgel into the edgemap using subpixel coordinates
  void insert(dbdet_edgel* e)
  {
    //determine appropriate cell to put this token into
    int xx = dbdet_round(e->pt.x());
    if(xx==0)
	xx++;
    if(xx>=width())
	xx=width()-1;
    int yy = dbdet_round(e->pt.y());
    if(yy==0)
	yy++;
    if(yy>=height())
	yy=height()-1;
    insert(e, xx, yy);
  }

  bool AlmostEqual(dbdet_edgemap & map, double angle_thresh, double strength_thresh) {

    bool ret = true;
    if (this==&map) {
      return true;
    }
    ret &= (this->num_edgels() == map.num_edgels());
    ret &= (this->width() == map.width());
    ret &= (this->height() == map.height());

    if (ret) {
      for (unsigned int xx=0; xx<width(); xx++) {
        for (unsigned int yy=0; yy<height(); yy++) {
          vcl_vector<dbdet_edgel*> cellA = this->cell(xx,yy);
          vcl_vector<dbdet_edgel*> cellB = map.cell(xx,yy);
          if (cellA.size() == cellB.size()) {
            for (unsigned int i=0; i<cellA.size(); i++) {
	      if(cellA[i]->pt != cellB[i]->pt) {
		return false;
	      }
	      if(vnl_math::abs(cellA[i]->tangent - cellB[i]->tangent) > angle_thresh) {
		return false;
	      }
              if(vnl_math::abs(cellA[i]->strength - cellB[i]->strength) > strength_thresh) {
		return false;
	      }
            }
          } else {
            return false;
          }
        }
      }
    } else {
	return false;
    }
    return true;
  }

  bool operator ==(dbdet_edgemap & map) {

    bool ret = true;
    if (this==&map) {
      return ret;
    }
    ret &= (this->num_edgels() == map.num_edgels());
    ret &= (this->width() == map.width());
    ret &= (this->height() == map.height());

    if (ret) {
      for (unsigned int xx=0; xx<width(); xx++) {
        for(unsigned int yy=0; yy<height(); yy++) {
          vcl_vector<dbdet_edgel*> cellA = this->cell(xx,yy);
          vcl_vector<dbdet_edgel*> cellB = map.cell(xx,yy);
          if (ret && cellA.size() == cellB.size()) {
            for (unsigned int i=0; i<cellA.size(); i++) {
              ret &= cellA[i]->pt == cellB[i]->pt;
              ret &= cellA[i]->tangent == cellB[i]->tangent;
              ret &= cellA[i]->strength == cellB[i]->strength;
              if (!ret)
                return false;
            }
          } else {
            return false;
          }
        }
      }
    }
    return ret;
  }
};


#endif // dbdet_edgemap_h
