// This is dbdet_appearance.h
#ifndef dbdet_appearance_h
#define dbdet_appearance_h
//:
//\file
//\brief The base class for storing edgel appearance info
//\author Amir Tamrakar
//\date 09/05/07
//
//\endverbatim

#include <vcl_string.h>
#include <vcl_iostream.h>
#include <vcl_sstream.h>

#include "dbdet_EMD.h"

//: Base class for image appearance on edgels
class dbdet_appearance
{
public:
  
  //: default constructor
  dbdet_appearance(){}

  //: copy constructor
  dbdet_appearance(const dbdet_appearance& /*other*/) {}

  //: destructor
  virtual ~dbdet_appearance(){}

  //: compute the distance between two appearance values
  virtual double dist(const dbdet_appearance& /*other*/)=0;

  virtual double value()=0;

  //: 
  virtual vcl_string print_info() const=0;

  //: Clone `this': creation of a new object and initialization. This is used so
  // that client classes like dbdet_edgel supports any type of dbdet_appearance.
  // See Prototype design pattern in Gamma et. al.
  virtual dbdet_appearance *clone() const = 0;
};


//: A derived class to store intensity as appearance
class dbdet_intensity: public dbdet_appearance
{
public:
  double val;      ///< the value of the intensity

  //: default constructor
  dbdet_intensity() : val(0.0){}
  
  //: constructor
  dbdet_intensity(double value) : val(value){}
  
  //: copy constructor
  dbdet_intensity(const dbdet_intensity & other) : dbdet_appearance(other) {val = other.val; }

  //: destructor
  ~dbdet_intensity(){}

  virtual double value(){ return val; }

  //: compute the distance between two appearance values (assumes that the operator '-' is defined for class A)
  virtual double dist(const dbdet_appearance& other)
  {
    return val - ((dbdet_intensity*)&other)->val;
  }

  //: return a string with the info of this appearance measure
  virtual vcl_string print_info() const { 
    vcl_stringstream ss;
    ss << val << '\0';
    return ss.str();
  }

  //: \see dbdet_appearance::clone
  virtual dbdet_appearance *clone() const {
    return new dbdet_intensity(*this);
  }
};

//: A derived class to store color as appearance
class dbdet_color: public dbdet_appearance
{
public:
  double c1, c2, c3;      ///< the value of the components

  //: default constructor
  dbdet_color() : c1(0.0), c2(0.0), c3(0.0){}
  
  //: constructor
  dbdet_color(double v1, double v2, double v3) : c1(v1), c2(v2), c3(v3){}
  
  //: copy constructor
  dbdet_color(const dbdet_color & other) : dbdet_appearance(other) {c1 = other.c1; c2 = other.c2; c3 = other.c3;}

  //: destructor
  ~dbdet_color(){}

  virtual double value(){ return c1; } //not meaningful

  //: compute the distance between two appearance values (assumes that the operator '-' is defined for class A)
  virtual double dist(const dbdet_appearance& other)
  {
    //euclidean distance
    dbdet_color* o = (dbdet_color*)&other;

    return vcl_sqrt((c1-o->c1)*(c1-o->c1) + (c2-o->c2)*(c2-o->c2) + (c3-o->c3)*(c3-o->c3));
  }

  //: return a string with the info of this appearance measure
  virtual vcl_string print_info() const{ 
    vcl_stringstream ss;
    ss << "(" << c1 << ", " << c2 << ", " << c3 << ")" << '\0';
    return ss.str();
  }

  //: \see dbdet_appearance::clone
  virtual dbdet_appearance *clone() const {
    return new dbdet_color(*this);
  }
};


//: A derived class to store gray dbdet_signature as appearance
class dbdet_gray_signature: public dbdet_appearance
{
public:
  dbdet_signature sig;      ///< the stored signature

  //: default constructor
  dbdet_gray_signature(){}
  
  //: constructor
  dbdet_gray_signature(dbdet_signature newsig) : sig(newsig){}
  
  //: copy constructor
  dbdet_gray_signature(const dbdet_gray_signature & other):dbdet_appearance(other) {sig = other.sig; }

  //: destructor
  ~dbdet_gray_signature(){}

  virtual double value(){ return 0.0; }

  //: compute the distance between two appearance values (assumes that the operator '-' is defined for class A)
  virtual double dist(const dbdet_appearance& other)
  {
    return vcl_fabs(sig - ((dbdet_gray_signature*)&other)->sig);
  }

  //: return a string with the info of this appearance measure
  virtual vcl_string print_info() const
  { 
    vcl_stringstream ss;
    ss << "[" ;
    for (int i=0; i<NBINS; i++)
      ss << sig.bins[i].weight << " ,";
      //ss << "(" << sig.bins[i].value << ", " << sig.bins[i].weight << ") ,";
    ss << "]";
    return ss.str();  
  }

  //: \see dbdet_appearance::clone
  virtual dbdet_appearance *clone() const {
    return new dbdet_gray_signature(*this);
  }
};

//: A derived class to store gray dbdet_signature as appearance
class dbdet_color_signature: public dbdet_appearance
{
public:
  dbdet_color_sig sig;      ///< the stored signature

  //: default constructor
  dbdet_color_signature(){}
  
  //: constructor
  dbdet_color_signature(dbdet_color_sig newsig) : sig(newsig){}
  
  //: copy constructor
  dbdet_color_signature(const dbdet_color_signature & other) : dbdet_appearance(other) {sig = other.sig; }

  //: destructor
  ~dbdet_color_signature(){}

  virtual double value(){ return 0.0; }

  //: compute the distance between two appearance values (assumes that the operator '-' is defined for class A)
  virtual double dist(const dbdet_appearance& other)
  {
    return vcl_fabs(sig - ((dbdet_color_signature*)&other)->sig);
  }

  //: return a string with the info of this appearance measure
  virtual vcl_string print_info() const
  { 
    vcl_stringstream ss;
    ss << "[" ;
    //for (int i=0; i<NBINS; i++)
    //  ss << sig.bins[i].weight << " ,";
      //ss << "(" << sig.bins[i].value << ", " << sig.bins[i].weight << ") ,";
    //ss << "]";
    return ss.str();  
  }

  //: \see dbdet_appearance::clone
  virtual dbdet_appearance *clone() const {
    return new dbdet_color_signature(*this);
  }
};

//: independent global function to compute distance between two appearances
inline double dbdet_appearance_dist(dbdet_appearance* app1, dbdet_appearance* app2)
{
  return app1->dist(*app2);
}

#endif // dbdet_appearance_h
