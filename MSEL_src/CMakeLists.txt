cmake_minimum_required(VERSION 2.8)
project( sel_linker )

FIND_PACKAGE(VXL)
#FIND_PACKAGE(VXD)
#FIND_PACKAGE(LEMS)


IF(VXL_FOUND)
	# include the standard VXL stuff
	#SET(VXL_PROVIDE_STANDARD_OPTIONS 1)
	INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)
	INCLUDE(${VXL_DIR}/VXLConfig.cmake)
	INCLUDE_DIRECTORIES( ${VXL_GEL_INCLUDE_DIR} )
	INCLUDE_DIRECTORIES( ${VXL_MUL_INCLUDE_DIR} )	
	INCLUDE_DIRECTORIES( ${VXL_BRL_INCLUDE_DIR} )	
	#MESSAGE("VXL_INCLUDE_DIR is ${VXL_CORE_INCLUDE_DIR}")
	MESSAGE("VXL_ROOT_SOURCE_DIR is ${VXL_ROOT_SOURCE_DIR}")

	#set( dborl_edge_det_link_sources
	#	dborl_edge_det_link_main.cxx
	#)

	set( MSEL_from_edge_sources
	  	MSEL_compute_contours_from_edges_main.cxx 
	)

	set( MSEL_from_img_sources
	  	MSEL_compute_contours_from_img_main.cxx 
	)

	#add_executable (dborl_compute_curve_frags ${dborl_edge_det_link_sources} )
	#target_link_libraries (dborl_compute_curve_frags MSEL_core) 

	add_executable (MSEL_compute_contours_from_edges ${MSEL_from_edge_sources} )
	target_link_libraries (MSEL_compute_contours_from_edges MSEL_core)

	add_executable (MSEL_compute_contours_from_img ${MSEL_from_img_sources} )
	target_link_libraries (MSEL_compute_contours_from_img MSEL_core)

	add_subdirectory( MSEL_core )
	
ENDIF(VXL_FOUND)
