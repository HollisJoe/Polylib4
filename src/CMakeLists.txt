
########### install files ###############

# install files 
if(SERIALTARGET)

 	include_directories(${PROJECT_SOURCE_DIR}/include)
 	include_directories(${TP_DIR}/include)
    if(USE_NPATCH_LIB)
 		include_directories(${NPT_DIR}/include)
 	endif()

    add_definitions("${REAL_OPT} ${NPT_OPT}")

	add_library(POLY Polylib.cxx 
 	     Polylib.cxx
     	 groups/PolygonGroup.cxx 
    	 groups/VTree.cxx 
     	 polygons/Triangle.cxx 
     	 polygons/NptTriangle.cxx 
     	 file_io/FileIO_func.cxx 
     	 file_io/PolygonIO.cxx 
     	 c_lang/CPolylib.cxx 
     	 f_lang/FPolylib.cxx 
     	 util/time.cxx
	)

	install(TARGETS POLY DESTINATION ${PL}/lib)

else()

 	include_directories(${PROJECT_SOURCE_DIR}/include)
 	include_directories(${MPI_DIR}/include)
 	include_directories(${TP_DIR}/include)
 	include_directories(${NPT_DIR}/include)

    add_definitions("${MPI_PL_OPT} ${REAL_OPT} ${NPT_OPT}")

	add_library(POLYmpi  Polylib.cxx
	 Polylib_MPI.cxx 
	 groups/PolygonGroup.cxx 
	 groups/PolygonGroup_MPI.cxx 
	 groups/VTree.cxx 
	 polygons/Triangle.cxx 
	 polygons/NptTriangle.cxx 
	 file_io/FileIO_func.cxx 
	 file_io/PolygonIO.cxx 
	 c_lang/CPolylib.cxx 
	 f_lang/FPolylib.cxx
	 util/time.cxx
	)

	install(TARGETS POLYmpi DESTINATION ${PL}/lib)
endif()


INSTALL(FILES ../include/Polylib_func.h ../include/Polylib.h
              ${PROJECT_BINARY_DIR}/include/polyVersion.h 
        DESTINATION ${PL}/include)

INSTALL(FILES ../include/c_lang/CPolylib.h
DESTINATION ${PL}/include/c_lang)

INSTALL(FILES ../include/f_lang/FPolylib.h
../include/f_lang/FPolylib_User.inc
../include/f_lang/FPolylib_User_f77.inc
../include/f_lang/FPolylib_define.inc
../include/f_lang/FPolylib_define_f77.inc
../include/f_lang/FPolylib_precision.inc
../include/f_lang/FPolylib_precision_def.inc
../include/f_lang/FPolylib_prototype.inc
../include/f_lang/FPolylib_prototype_f77_double.inc
../include/f_lang/FPolylib_prototype_f77_float.inc
../include/f_lang/FPolylib_struct.inc
DESTINATION ${PL}/include/f_lang)

INSTALL(FILES ../include/common/BBox.h
../include/common/PolylibCommon.h ../include/common/PolylibDefine.h
../include/common/PolylibStat.h ../include/common/tt.h
../include/common/Vec2.h ../include/common/Vec3.h DESTINATION ${PL}/include/common)

INSTALL(FILES ../include/file_io/FileIO_func.h ../include/file_io/PolygonIO.h
DESTINATION ${PL}/include/file_io)

INSTALL(FILES ../include/groups/PolygonGroup.h
../include/groups/VTree.h DESTINATION ${PL}/include/groups)

INSTALL(FILES ../include/polygons/NptTriangle.h
../include/polygons/Triangle.h DESTINATION ${PL}/include/polygons)

 
INSTALL(FILES ../include/util/time.h DESTINATION ${PL}/include/util)


SET(EXTRA_DIST "${EXTRA_DIST} Makefile_hane depend.inc")

###
