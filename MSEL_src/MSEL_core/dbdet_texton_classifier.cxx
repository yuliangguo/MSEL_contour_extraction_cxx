#include <vcl_cassert.h>
#include <vcl_cstring.h>
#include <vcl_limits.h>
#include <vnl/vnl_file_matrix.h>
#include "dbdet_texton_classifier.h"
#include "dbdet_filter_util.h"

dbdet_texton_classifier::dbdet_texton_classifier(const char * filename)
{
  classes = static_cast<vnl_matrix <double> > (vnl_file_matrix<double>(filename));
}

int dbdet_texton_classifier::numClasses()
{
  return classes.cols();
}

vnl_matrix<unsigned> dbdet_texton_classifier::classify(vcl_vector<vil_image_view<double> > decomposed)
{
  assert(decomposed.size() == classes.rows() && decomposed.size() > 0);

  unsigned xs = decomposed[0].ni() * decomposed[0].nj();
  unsigned ys = classes.cols();

  vnl_matrix<double> dec(decomposed.size(), xs);
  for (unsigned i = 0; i < decomposed.size(); ++i)
    dec.set_row(i, (decomposed[i]).top_left_ptr());

  vcl_vector<double> x2(xs, 0.);
  vcl_vector<double> y2(ys, 0.);

  for (unsigned i = 0; i < dec.rows(); ++i)
    for (unsigned j = 0; j < dec.cols(); ++j)
      x2[j] += dec(i, j) * dec(i, j);

  for (unsigned i = 0; i < classes.rows(); ++i)
    for (unsigned j = 0; j < classes.cols(); ++j)
      y2[j] += classes(i, j) * classes(i, j);

  vnl_matrix<unsigned> classified(decomposed[0].ni(), decomposed[0].nj());
  for (unsigned i = 0; i < xs; ++i)
  {
    vcl_vector<double> d(ys, 0);
    unsigned idxj = i / classified.rows();
    unsigned idxi = i % classified.rows();


    for (unsigned k = 0; k < dec.rows(); ++k)
      for (unsigned j = 0; j < ys; ++j)
        d[j] += dec(k, i) * classes(k, j);

    double dist = vcl_numeric_limits<double>::max();
    for (unsigned j = 0; j < ys; ++j)
    {
      double tmp = -2. * d[j] + x2[i] + y2[j];
      if (tmp < dist)
      {
        dist = tmp;
        
        classified(idxi, idxj) = j;
      }
    }
  }
  return classified;
}
