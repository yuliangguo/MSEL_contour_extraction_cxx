// This is dbgl_diffgeom.h
#ifndef dbgl_diffgeom_h_
#define dbgl_diffgeom_h_
//:
//\file
//\brief Differential geometry estimation for curves (extensions to VXD)
//\author Based on original code by  Ricardo Fabbri (rfabbri), Brown University  (rfabbri.github.io)
//\date Wed May 25 15:53:23 BRT 2016
//

#include "dbdet_edgel.h"
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>

//: Same as the vxd/bgld version, but for dbdet_edgel_chain
void 
dbgl_compute_curvature(
    const dbdet_edgel_chain &c, 
    vnl_vector<double> *k
    );

//: Same as the vxd/bgld version, but for dbdet_edgel_chain
void 
dbgl_compute_normals(
    const dbdet_edgel_chain &c, 
    vcl_vector< vnl_vector_fixed<double, 2> > *n
    );

#endif // dbgl_diffgeom_h

