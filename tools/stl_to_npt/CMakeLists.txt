
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${MPI_DIR}/include)
include_directories(${TP_DIR}/include)

add_definitions("${MPI_PL_OPT} ${REAL_OPT} ${NPT_OPT}")


link_directories(${PROJECT_BINARY_DIR}/src/${CMAKE_CONFIGURATION_TYPES})
link_libraries(${PL_LIB})

if(USE_NPATCH_LIB)
    include_directories(${NPT_DIR}/include)
    link_directories(${NPT_DIR}/lib)
    link_libraries(${NPT_LIB})
endif()

link_directories(${TP_DIR}/lib)
link_libraries(${TP_LIB})

if(enable_mpi STREQUAL msmpi)
    link_directories(${MPI_DIR}/lib/x64)
    link_libraries(${MPI_LIB})
endif()


if(SERIALTARGET)
set(stl_to_npt_SRCS
   main.cxx
)
else()
set(stl_to_npt_SRCS
   main.cxx
   gather_polygons_to_npt.cxx
   StlToNpt.h
)
endif()

add_executable(stl_to_npt ${stl_to_npt_SRCS})

install(TARGETS stl_to_npt DESTINATION ${PL}/bin)
