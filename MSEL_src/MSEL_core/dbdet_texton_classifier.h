#include <vnl/vnl_matrix.h>
#include <vcl_vector.h>
#include <vil/vil_image_view.h>

#ifndef _DBDET_TEXTON_CLASSIFIER_H_
#define _DBDET_TEXTON_CLASSIFIER_H_

class dbdet_texton_classifier {
public:
  vnl_matrix<double> classes;


  dbdet_texton_classifier(const char * filename);

  int numClasses();

  vnl_matrix<unsigned> classify(vcl_vector<vil_image_view<double> > decomposed);
};

#endif
