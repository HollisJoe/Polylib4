
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include/f_lang)
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


set(f_interface_SRCS
   main.F90
)

add_executable(f_interface ${f_interface_SRCS})

# 特殊な指定が必要
#   ${ADD_TEST_LINK_LIB}: 複数のライブラリがある場合、先頭の-lは省略、２番め以降は-lをつける   
target_link_libraries(f_interface  ${ADD_TEST_LINK_LIB})

#--- 上記でlinkまでを行う -------
#   テストは手動で行う
#   make test 未対応