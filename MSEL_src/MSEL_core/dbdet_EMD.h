// This is brcv/seg/dbdet/algo/dbdet_EMD.h
#ifndef dbdet_EMD_h
#define dbdet_EMD_h
//:
//\file
//\brief An implementation of the EMD algorithm to compute distance between two signatures
//\author Amir Tamrakar
//\date 08/26/07
//
//\verbatim
//  Modifications
//\endverbatim

#include <vcl_cmath.h>
#include <vcl_algorithm.h>

//#include "dbdet_BS.h"

//for grayscale
#define NBINS 8
#define GAMMA 14.0
#define MAXCLUSTERS 30

typedef struct bin_struct{
  double wsum;
  double weight;
  double value;
} dbdet_bin;

//for color
typedef struct
{
  int from;             /* Feature number in signature 1 */
  int to;               /* Feature number in signature 2 */
  double amount;        /* Amount of flow from "from" to "to" */
} dbdet_flow;

//forward declaration
class dbdet_signature;
class dbdet_color_sig;

//------------------------------------------------------------------------------------------
// Grayscale histogram distance functions
double dbdet_gray_EMD(const dbdet_bin dirt[], const dbdet_bin hole[]);
double dbdet_chi_sq_dist(const dbdet_bin hist1[], const dbdet_bin hist2[]);
double dbdet_bhat_dist(const dbdet_bin hist1[], const dbdet_bin hist2[]);
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------
//Color Histogram distance functions
double dbdet_color_EMD(dbdet_color_sig *sig1, dbdet_color_sig *sig2, dbdet_flow* flow=NULL, int* flowsize=NULL);
double dbdet_color_chi_sq_dist(const dbdet_color_sig &sig1, const dbdet_color_sig &sig2);
double dbdet_color_bhat_dist(const dbdet_color_sig &sig1, const dbdet_color_sig &sig2);
//------------------------------------------------------------------------------------------

//: This code implements the basic Binary Splitting algorithm described in 
// the 1991 IEEE Trans. on Sig. Proc. article "Color Quantization of Images"
// by Michael Orchard and Charles Bouman, pp. 2677-90.
//
// The input is a 1-D array of data points.  The #define'd constant CDIM 
// holds the number of dimensions.  The output is a set of cluster centers
// and is also a 1-D array.  The number of clusters is also returned.
void bs(float *, int, int, int *, float **, int *);

//: This class represents an intensity distribution as a signature (an adaptive histogram)
//  The bin centers of this adaptive histogram are dynamically determined. The weight parameter stores
//  mass of that bin.
class dbdet_signature
{
public:
  dbdet_bin bins[NBINS]; //fixed number of bins for now (this is not strictly required)

  dbdet_signature()
  {
    for (int i=0; i<NBINS; i++){
      bins[i].value=0; 
      bins[i].weight=0;
      bins[i].wsum=0;
    }
  }
  ~dbdet_signature(){}

  //: The EMD dist is the default distance between two signatures
  double operator-(const dbdet_signature& sig) const 
  { 
    //return dbdet_chi_sq_dist(bins, sig.bins); 
    return dbdet_gray_EMD(bins, sig.bins); 
  }

  //: This operator sums two signatures
  void operator +=(dbdet_signature const& sig) 
  { 
    for (int i=0; i<NBINS; i++){
      bins[i].weight += sig.bins[i].weight;
      bins[i].weight /= 2.0;
      bins[i].wsum += sig.bins[i].wsum;
      //renormalize
      bins[i].value = sig.bins[i].wsum/bins[i].weight; 
    }
  }

};

//: This class represents the Lab color distribution as a signature (an adaptive histogram)
//  The bin centers of this adaptive histogram are dynamically determined. The weight parameter stores
//  mass of that bin.
class dbdet_color_sig
{
public:
  int n; //number of features
  double Features[MAXCLUSTERS]; //fixed number of features for now (this is not strictly required)
  double Weights[MAXCLUSTERS]; 

  dbdet_color_sig(){}
  ~dbdet_color_sig(){}

  //: The EMD dist is the default distance between two signatures
  double operator-(const dbdet_color_sig& /*sig1*/) const 
  { 
    //return dbdet_chi_sq_dist(bins, sig1.bins); 
    //return dbdet_color_EMD(*this, sig1, NULL, NULL); 
    return 0;
  }

};


#endif // dbdet_EMD_h
