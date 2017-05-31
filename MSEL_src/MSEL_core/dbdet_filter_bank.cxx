#include <vcl_iostream.h>
#include <vcl_fstream.h>
#include <vcl_sstream.h>
#include <vcl_cassert.h>
#include <vnl/vnl_file_matrix.h>
#include <vil/vil_convert.h>

#include "dbdet_filter_util.h"
#include "dbdet_filter_bank.h"

dbdet_filter_2d::dbdet_filter_2d(const char * fname) {
  m = static_cast<vnl_matrix<double> > (vnl_file_matrix<double>(fname));
  assert(m.cols() % 2 == 1 && m.rows() % 2 == 1);
}

bool dbdet_filter_2d::isEmpty()
{
  return (m.cols() == 0 && m.rows() == 0);
}

int dbdet_filter_2d::size()
{
  return m.cols();
}

double dbdet_filter_2d::flipped(int i, int j)
{
  return m(m.rows() - i - 1, m.cols() - j - 1);
}

vil_image_view<double> dbdet_filter_2d::applyPadded(vil_image_view<vxl_byte> image, int border)
{
  vil_image_view<double> ret = vil_image_view<double>(image.ni() - 2 * border, image.nj() - 2 * border);

  int halfC = m.cols() / 2;
  int halfR = m.rows() / 2;

  for (int i = 0; i < ret.ni(); ++i)
  {
    int ib = i + border;
    for (int j = 0; j < ret.nj(); ++j)
    {
      int jb = j + border;
      ret(i, j) = 0.0;
      for (int ii = -halfR; ii <= halfR; ++ii)
        for (int jj = -halfC; jj <= halfC; ++jj)
          ret(i, j) += image(ib + ii, jb + jj) * flipped(ii + halfR, jj + halfC);
    }
  }
  return ret;
}

vil_image_view<double> dbdet_filter_2d::applyPadded13(vil_image_view<vxl_byte> image, int border)
{
  vil_image_view<double> ret = vil_image_view<double>(image.ni() - 2 * border, image.nj() - 2 * border);

  const int halfC = 6;
  const int halfR = 6;

  for (int i = 0; i < ret.ni(); ++i)
  {
    int ib = i + border;
    for (int j = 0; j < ret.nj(); ++j)
    {
      int jb = j + border;
      ret(i, j) = 0.0;
      for (int ii = -halfR; ii <= halfR; ++ii)
        for (int jj = -halfC; jj <= halfC; ++jj)
          ret(i, j) += image(ib + ii, jb + jj) * flipped(ii + halfR, jj + halfC);
    }
  }
  return ret;
}

vil_image_view<double> dbdet_filter_2d::applyPadded19(vil_image_view<vxl_byte> image, int border)
{
  vil_image_view<double> ret = vil_image_view<double>(image.ni() - 2 * border, image.nj() - 2 * border);

  const int halfC = 9;
  const int halfR = 9;

  for (int i = 0; i < ret.ni(); ++i)
  {
    int ib = i + border;
    for (int j = 0; j < ret.nj(); ++j)
    {
      int jb = j + border;
      ret(i, j) = 0.0;
      for (int ii = -halfR; ii <= halfR; ++ii)
        for (int jj = -halfC; jj <= halfC; ++jj)
          ret(i, j) += image(ib + ii, jb + jj) * flipped(ii + halfR, jj + halfC);
    }
  }
  return ret;
}




dbdet_filter_bank::dbdet_filter_bank(vcl_string baseDir)
{
  filtersMaxSize = 0;
  vcl_ifstream filterFiles((baseDir + "/filter_list.txt").c_str());
  if (filterFiles.good())
  {
    vcl_string line;
    while (!filterFiles.eof())
    {
      vcl_getline(filterFiles, line);
      if(line.length() > 0) 
      {
        line = baseDir + "/" + line;
        dbdet_filter_2d f = dbdet_filter_2d(line.c_str());
        assert(!f.isEmpty());
        filters.push_back(f);
        filtersMaxSize = f.size() > filtersMaxSize ? f.size() : filtersMaxSize;
      }
    }
  }
}

int dbdet_filter_bank::numFilters()
{
  return filters.size();
}

vcl_vector<vil_image_view<double> > dbdet_filter_bank::decompose(vil_image_view<vil_rgb<vxl_byte> > image)
{
  vcl_vector<vil_image_view<double> > ret;
  ret.resize(filters.size());
  int padSize = filtersMaxSize / 2;
  vil_image_view<vxl_byte> grey_img;
  //Using luma weights to match matlab
  vil_convert_rgb_to_grey(image, grey_img, 0.2989, 0.5870, 0.1140);
  grey_img = padReflect(grey_img, padSize);

  for (int i = 0; i < filters.size(); ++i)
  {
    switch (filters[i].size())
    {
    case 13:
      ret[i] = filters[i].applyPadded13(grey_img, padSize);
      break;
    case 19:
      ret[i] = filters[i].applyPadded19(grey_img, padSize);
      break;
    default:
      ret[i] = filters[i].applyPadded(grey_img, padSize);
    }
  }
  return ret;
}
