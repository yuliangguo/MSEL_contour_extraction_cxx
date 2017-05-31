# MSEL_contour_extraction_cxx

Author of this Release Package: 
	Yuliang Guo (yuliang_guo@brown.edu)
	This is multi-stage approach in extracting curve fragments features from image.
	This package includes research code still under development. There are a lot redundent code included.
	The evaluation result is not the same reported in the published papers.

Reference: 
	"A Multi-Stage Approach to Curve Extraction", Y.Guo, N.Kumar, M.Narayanan and B.Kimia, ECCV 2014
	"On Evaluating Methods for Recovering Image Curve Fragments", Y.Guo, B.Kimia, CVPRW 2012
    	"No grouping left behind: From edges to curve fragments, Tamrakar and Kimia, ICCV 2007"


### 1. Download VXL

```bash
  mkdir vxl
  git clone https://github.com/vxl/vxl.git vxl
```

### 2. Compile VXL
```bash
  mkdir vxl-bin
  cd ./vxl-bin
  ccmake -D CMAKE_BUILD_TYPE=Release -D BUILD_SHARED_LIBS=ON -D VNL_CONFIG_LEGACY_METHODS=ON ../vxl
  # press 'c' (configure) multiple times until 'g' (generate) appears
  make -j4 -k   # compile in parallel and keep going past errors

  # Don't worry about errors at this point. We will not use everything.
```

### 3. compile MSEL
```bash
	mkdir MSEL_bin
	cd ./MSEL_bin
	ccmake -D CMAKE_BUILD_TYPE=Release ../MSEL_src
	# type in the path of VXL_DIR as the path of vxl-bin folder made in step 2.
  	# press 'c' (configure) multiple times until 'g' (generate) appears
	make -j4 -k
	# there are two excutables generated in MSEL_bin
```
### 4. Use of MSEL contour extraction

	# Extraction contours from image 
		e.g.	./MESL_bin/MSEL_compute_contours_from_img 10081.jpg 10081.cem 200 1.5 1
	# This is a combination of edge-detection and contour extraction
	# the edge detection integrated referst to: "No grouping left behind: From edges to curve fragments, Tamrakar and Kimia, ICCV 2007"
	# Usage: ./MESL_bin/MSEL_compute_contours_from_img input_img_file output_cem_file nContours edge_sigma edge_thresh
	# 	input_img_file: path of the input image file. Input image must be color.
	#	output_cem_file: path of the output contour file, we define this type of file as ".cem".
	#	nContours: number of contours to be outputed after ranking. If this is set to 0 or kept blank, it will output all the generated contours.
	#	edge_sigma:  sigma parameter for the edge detection, deciding the scale of edge.
	#	edge_thresh: gradient thresh for the edge detection, deciding the scale of edge.

	# Extraction contours from edges 
		e.g.	./MESL_bin/MSEL_compute_contours_from_edges 10081.jpg 10081.edg 10081.cem 200
	# Suppose edges can be detected from third-part softwares, and can be saved in the same format as "10081.edg"
	# Usage: ./MESL_bin/MSEL_compute_contours_from_img input_img_file input_edg_file output_cem_file nContours
	# 	input_img_file:	path of the input image file. Input image must be color.
	# 	input_edg_file: path of the input edgemap file, which we define the type as ".edg".
	#	output_cem_file: path of the output contour file, we define this type of file as ".cem".
	#	nContours: number of contours to be outputed after ranking. If this is set to 0 or kept blank, it will output all the generated contours.

### 5. Visualize contours using Matlab

 	use Matlab, demo_vis_io.m

