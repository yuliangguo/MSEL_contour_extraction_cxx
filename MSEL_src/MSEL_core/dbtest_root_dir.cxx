// This is basic/dbtest/dbtest_root_dir.cxx


//:
// \file



#include "dbtest_root_dir.h"

//#include <vcl_cstdlib.h>
#include <vcl_iostream.h>

// The following should have been created automatically by the
// configuration scripts from dbtest_where_root_dir.h.in
// We need to check for its existence and if it doesn't exist - do something else.

#ifdef DBTEST_WHERE_ROOT_DIR_H_EXISTS
#include "dbtest_where_root_dir.h"

//: Return source root directory of LEMSVXL source
vcl_string dbtest_root_dir()
{
  return vcl_string(DBTEST_SOURCE_ROOT_DIR);
}
#else
//: Return source root directory of LEMSVXL source
vcl_string dbtest_root_dir()
{
  vcl_cerr << "ERROR: dbtest_root_dir() Unable to retrieve root directory.\n"
           << "Missing header file dbtest_where_root_dir.h. Please check configuration of dbtest.\n";
  return vcl_string("");
}

#endif
