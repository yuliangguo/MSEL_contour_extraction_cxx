// This is basic/dbtest/dbtest_root_dir.h

#ifndef dbtest_root_dir_h_
#define dbtest_root_dir_h_

//:
// \file
// \brief Function to return root directory of lemsvxl source.
//        This is similar to testlib_root_dir.h in VXL
// \author Nhon Trinh (ntrinh@lems.brown.edu)
// \date Oct 30, 2009
// \verbatim
// \endverbatim

#include <vcl_string.h>

//: Return source root directory of LEMSVXL
//  If the file dbtest_where_root_dir.h has been automatically generated
//  during configuration (which will happen with cmake) then the
//  appropriate source directory will be returned.
//
//  If another build system is used in which this is not created,
//  the function will return empty space, or "".


vcl_string dbtest_root_dir();



#endif // testlib_root_dir_h_
