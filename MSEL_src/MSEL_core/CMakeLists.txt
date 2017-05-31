
	#INCLUDE_DIRECTORIES( ${VXL_GEL_INCLUDE_DIR} )
	#INCLUDE_DIRECTORIES( ${VXL_MUL_INCLUDE_DIR} )
	#MESSAGE("VXL_INCLUDE_DIR is ${VXL_CORE_INCLUDE_DIR}")
	#MESSAGE("VXL_GEL_INCLUDE_DIR is ${VXL_GEL_INCLUDE_DIR}")
	#MESSAGE("VXL_MUL_INCLUDE_DIR is ${VXL_MUL_INCLUDE_DIR}")

	configure_file(
	  ${CMAKE_SOURCE_DIR}/MSEL_core/dbtest_where_root_dir.h.in
	  ${CMAKE_BINARY_DIR}/MSEL_core/dbtest_where_root_dir.h ESCAPE_QUOTES
	)

	find_path(DBTEST_WHERE_ROOT_DIR_H_DIRECTORY
		  dbtest_where_root_dir.h
		  ${CMAKE_BINARY_DIR}/MSEL_core
		  DOC "Path to directory containing the file dbtest_where_root_dir.h")

	if (DBTEST_WHERE_ROOT_DIR_H_DIRECTORY)
	  add_definitions(-DDBTEST_WHERE_ROOT_DIR_H_EXISTS)
	  include_directories(${DBTEST_WHERE_ROOT_DIR_H_DIRECTORY})
	endif(DBTEST_WHERE_ROOT_DIR_H_DIRECTORY)


	set( MSEL_core_sources



		# From bgld and bnld in vxd
		bgld_diffgeom.h 			bgld_diffgeom.cxx
		bgld_eulerspiral_base.h
		bgld_eulerspiral.h			bgld_eulerspiral.cxx
		bgld_arc.h 				bgld_arc.cxx
		bgld_arc_base.h 			bgld_arc_base.cxx
		bgld_arc_algo.h 			bgld_arc_algo.cxx
		bgld_biarc.h 				bgld_biarc.cxx
		bgld_curve_smoothing.h			bgld_curve_smoothing.cxx
		bgld_param_curve.h
		bnld_fresnel.h				bnld_fresnel.cxx

		# From breye1d in vxd
		bpro1_parameters.h			bpro1_parameters.cxx			bpro1_parameters_sptr.h
		bpro1_parameters.hxx		
		bpro1_process.h 			bpro1_process.cxx			bpro1_process_sptr.h
		bpro1_storage.h			 	bpro1_storage.cxx			bpro1_storage_sptr.h
		vidpro1_vsol2D_storage.h		vidpro1_vsol2D_storage.cxx		vidpro1_vsol2D_storage_sptr.h
		vidpro1_image_storage.h			vidpro1_image_storage.cxx		vidpro1_image_storage_sptr.h
		vidpro1_save_cem_process.h		vidpro1_save_cem_process.cxx
		bsold_file_io.h				bsold_file_io.cxx

		# Templates from breye1d in vxd
		bpro1_parameters+bool-.cxx
		bpro1_parameters+bpro_filepath-.cxx
		bpro1_parameters+double-.cxx
		bpro1_parameters+vcl_string-.cxx
		bpro1_parameters+unsigned-.cxx
		bpro1_parameters+int-.cxx
		vbl_smart_ptr+bpro1_storage-.cxx
		vbl_smart_ptr+bpro1_process-.cxx
		vbl_smart_ptr+bpro1_parameters-.cxx				
		vsl_map_io+vcl_string.vcl_vector+vsol_spatial_object_2d_sptr--.cxx
		vsl_map_io+vcl_string.vcl_vector+double--.cxx
		vsl_binary_loader+bpro1_storage-.cxx



		# dbdet codes from lemsvxl
		dbgl_diffgeom.h				dbgl_diffgeom.cxx	
		dbtest_root_dir.h			dbtest_root_dir.cxx
		dbdet_appearance.h
		dbdet_cem_file_io.h			dbdet_cem_file_io.cxx
		dbdet_CFTG.h
		dbdet_curve_fragment_graph.h
		dbdet_curvelet.h			dbdet_curvelet.cxx
		dbdet_curvelet_map.h			dbdet_curvelet_map.cxx
		dbdet_curve_model.h			dbdet_curve_model.cxx
		dbdet_edgel.h				dbdet_edgel.cxx
		dbdet_edgel_link_graph.h
		dbdet_edgemap.h				dbdet_edgemap.cxx			dbdet_edgemap_sptr.h
		dbdet_edgemap_storage.h			dbdet_edgemap_storage.cxx		dbdet_edgemap_storage_sptr.h
		dbdet_EHT.h
		dbdet_EMD.h				dbdet_EMD.cxx
		dbdet_hyp_tree_graph.h
		dbdet_load_edg.h			dbdet_load_edg.cxx
		dbdet_load_edg_process.h		dbdet_load_edg_process.cxx
		dbdet_postprocess_contours.h		dbdet_postprocess_contours.cxx
		dbdet_save_cem_process.h		dbdet_save_cem_process.cxx
		dbdet_sel_base_CFTG_algo.cxx
		dbdet_sel_base.h			dbdet_sel_base.cxx

		dbdet_sel.cxx
		dbdet_sel+dbdet_CC_curve_model_3d-.cxx
		dbdet_sel+dbdet_CC_curve_model-.cxx
		dbdet_sel+dbdet_CC_curve_model_new-.cxx
		dbdet_sel+dbdet_CC_curve_model_perturbed-.cxx
		dbdet_sel+dbdet_ES_curve_model-.cxx
		dbdet_sel+dbdet_ES_curve_model_perturbed-.cxx
		dbdet_sel+dbdet_linear_curve_model-.cxx
		dbdet_sel+dbdet_simple_linear_curve_model-.cxx
		dbdet_sel.h									dbdet_sel_sptr.h
		dbdet_sel_process.h			dbdet_sel_process.cxx
		dbdet_sel_storage.h			dbdet_sel_storage.cxx			dbdet_sel_storage_sptr.h
		dbdet_sel_utils.h


		dbdet_contour_breaker.h			dbdet_contour_breaker.cxx
		dbdet_contour_ranker_process.h		dbdet_contour_ranker_process.cxx
		dbdet_curve_fragment_cues.h		dbdet_curve_fragment_cues.cxx
		dbdet_curve_fragment_ranker.h
		dbdet_edge_appearance_util.h		dbdet_edge_appearance_util.cxx
		dbdet_filter_bank.h			dbdet_filter_bank.cxx
		dbdet_filter_util.h
		dbdet_texton_classifier.h		dbdet_texton_classifier.cxx
		dbdet_yuliang_features.h		dbdet_yuliang_features.cxx
		dbdet_contour_breaker_geometric_process.h 	dbdet_contour_breaker_geometric_process.cxx
		dbdet_graphical_model_contour_merge_process.h	dbdet_graphical_model_contour_merge_process.cxx
		dbdet_graphical_model_contour_merge.h		dbdet_graphical_model_contour_merge.cxx

		dbdet_gaussian_kernel.h
		dbdet_interp_kernel.h
		dbdet_kernel.h
		dbdet_nms.h				dbdet_nms.cxx
		dbdet_subpix_convolution.h
		dbdet_third_order_color_detector.h	dbdet_third_order_color_detector.cxx
		dbdet_third_order_color_edge_detector_process.h	dbdet_third_order_color_edge_detector_process.cxx


		# Templates from dbdet in lemsvxl
		vbl_smart_ptr+dbdet_edgemap-.cxx
		vbl_smart_ptr+dbdet_sel_base-.cxx

	)
	add_library (MSEL_core ${MSEL_core_sources} )
	target_link_libraries (MSEL_core pdf1d vdgl vil_algo bil_algo brip) 


