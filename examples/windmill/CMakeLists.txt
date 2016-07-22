
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${MPI_DIR}/include)
include_directories(${TP_DIR}/include)
include_directories(./)

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


set(windmill_SRCS
   main.cxx
)

add_executable(windmill ${windmill_SRCS})
