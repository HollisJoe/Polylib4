###############################################################################
##
## Polylib - Polygon Management Library
##
## Copyright (c) 2010-2011 VCAD System Research Program, RIKEN.
## All rights reserved.
##
## Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
## All rights reserved.
##
###############################################################################

cmake_minimum_required(VERSION 3.0)

#enable_testing()

enable_language(CXX)
if(enable_test_c STREQUAL yes)
	set(TEST_C TRUE)
        enable_language(C)
endif()
if(enable_test_fortran STREQUAL yes)
	set(TEST_FORTRAN TRUE)
        enable_language(Fortran)
endif()

project(Polylib)
# Polylib install dir.
message(CMAKE_INSTALL_PREFIX: "${CMAKE_INSTALL_PREFIX}")
set(PL "${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}")

message(PROJECT_NAME      : "${PROJECT_NAME}")
message(PROJECT_SOURCE_DIR: "${PROJECT_SOURCE_DIR}")
message(PROJECT_BINARY_DIR: "${PROJECT_BINARY_DIR}")
message(CMAKE_INSTALL_PREFIX: "${CMAKE_INSTALL_PREFIX}")
message("Polylib install dir : ${PL}")

# Version info.
set(VERSION "4.0.0")
set( ${PROJECT_NAME}_MAJOR_VERSION 4)
set( ${PROJECT_NAME}_MINOR_VERSION 0)
set( ${PROJECT_NAME}_PATCH_LEVEL 0)
set(PACKAGE_VERSION ${VERSION})


# Revision No.
set(PL_REVISION "20160706_1200")

# Compile settting
#    CMAKE_CXX_COMPILER is not available on Windows. 
#    if you set CMAKE_CXX_COMPILER, but it is ignored.
#    CXX is reserved word.  it is used at enable_language().
#    so, use PL_CXX.  

set(CMAKE_CXX_COMPILER ${PL_CXX})
message(PL_CXX: "${PL_CXX}")
message(CMAKE_CXX_COMPILER: "${CMAKE_CXX_COMPILER}")

# for tests/c_interface or f_interface
set(ADD_TEST_LINK_LIB ${ADD_TEST_LINK_LIB})
message(ADD_TEST_LINK_LIB: "${ADD_TEST_LINK_LIB}")

#
# if a wrapper compiler is used, "enable_mpi" is empty, but parallel.
# wrappers like mpicxx, mpicc, mpiCC, mpic++, mpiFCCpx are taking into
# consideration, if not add in case.
# CC is for Hopper

if(${PL_CXX} STREQUAL mpicxx)
	set(enable_mpi "mpi")
	message("CMAKE_CXX_COMPILER=mpicxx")

elseif(${PL_CXX} STREQUAL mpicc)
	set(enable_mpi "mpi")
	message("CMAKE_CXX_COMPILER=mpicc")

elseif(${PL_CXX} STREQUAL mpiCC)
	set(enable_mpi "mpi")
	message("CMAKE_CXX_COMPILER=mpiCC")

elseif(${PL_CXX} STREQUAL mpic++)
	set(enable_mpi "mpi")
	message("CMAKE_CXX_COMPILER=mpic++")

elseif(${PL_CXX} STREQUAL mpiFCCpx)
	set(enable_mpi "mpi")
	message("CMAKE_CXX_COMPILER=mpiFCCpx")

elseif(${PL_CXX} STREQUAL CC)
	set(enable_mpi "mpi")
	message("CMAKE_CXX_COMPILER=CC")
else()
	set(enable_mpi "none")
	message("CMAKE_CXX_COMPILER is others")
endif()

# mpich //'-Dmpich=set_openMPI' '-Dwith_mpich=dir' option.
# openMPI //'-Dopmi=set_openMPI' '-Dwith_ompi=dir' option.
if(with_mpich)
	set(MPI_DIR ${with_mpich})
	set(enable_mpi "mpich")
elseif(with_ompi)
	set(MPI_DIR ${with_ompi})
	set(enable_mpi "ompi")
elseif(with_msmpi)
	set(MPI_DIR ${with_msmpi})
	set(enable_mpi "msmpi")
endif()

#warning
if(NOT enable_mpi STREQUAL mpi)
	message("Warning : Neither MPICH nor OpenMPI nor MSMPI library is not specified.")
	message("If compilation error occurs about MPI,specify -Dwith_mpich, -Dwith_ompi or -Dwith_msmpi option.")
endif()

# test execution command

# mpi execution command & mpi compile option
# 	${MPI_EXEC} ${MPI_PREFIX_NUM_PROC} ${MPI_TEST_NUM_PROC}
#	add -DMPI_PL for safe

if(enable_mpi STREQUAL mpich)
	set(MPI_PL_OPT "-DMPI_PL")
	set(TEST_RUN_MPI "test_run_mpiexec.sh")
#message("mpi success1")
elseif(enable_mpi STREQUAL ompi)
	set(MPI_PL_OPT "-DMPI_PL")
	set(TEST_RUN_MPI "test_run_mpirun.sh")
elseif(enable_mpi STREQUAL msmpi)
	set(MPI_PL_OPT "-DMPI_PL")
	set(MPI_LIB "msmpi")
	set(TEST_RUN_MPI "test_run_mpiexec.sh")
#message("mpi success2")
elseif(enable_mpi STREQUAL mpi)
	set(MPI_PL_OPT "-DMPI_PL")
	if(PL_CXX STREQUAL mpiFCCpx)
		set(TEST_RUN_MPI "test_run_mpiexec.sh")
	else()
		message("Warning: TEST_RUN_MPI is not specified.")
		# message("Tests will be error.")
	endif()
endif()

message(enable_mpi: "${enable_mpi}")
message(MPI_PL_OPT: "${MPI_PL_OPT}")

# TextParser environment

if(with_parser)
	set(TP_DIR ${with_parser})
else()
	message("Error: TextParser library is not specified.")
	message("specify -Dwith_parser option.")
endif()

# Npatch library

if(with_npatch)
	set(USE_NPATCH_LIB TRUE)
	set(NPT_OPT "-DUSE_NPATCH_LIB")
	set(NPT_DIR ${with_npatch})
	set(NPT_LIB "Npatch")
else()
#	message("Error: Nagata patch library is not specified.")
#	message("Specify -Dwith_npatch option.")
	set(NPT_OPT "")
	set(NPT_DIR "" )
	set(NPT_LIB "")
endif()
message(NPT_OPT: "${NPT_OPT}")
message(NPT_DIR: "${NPT_DIR}")
message(NPT_LIB: "${NPT_LIB}")


if(enable_mpi STREQUAL none)
	set(SERIALTARGET TRUE)
	set(PL_LIB "POLY")
	set(TP_LIB "TP")
else()
	set(SERIALTARGET FALSE)
	set(PL_LIB "POLYmpi")
	set(TP_LIB "TPmpi")
endif()

# Type of REAL

if(with_real STREQUAL double)
	add_definitions(-D_REAL_IS_DOUBLE_)
	set(REAL_OPT "${REAL_OPT} -D_REAL_IS_DOUBLE_")
endif()	


##Checks for header files.
#-----------------------------------------------------------------
include(CheckFunctionExists)
include(CheckIncludeFiles)

check_function_exists(gettimeofday HAVE_GETTIMEOFDAY) 
check_include_files(inttypes.h HAVE_INTTYPES_H)
check_include_files(memory.h HAVE_MEMORY_H)
check_function_exists(memset HAVE_MEMSET)
check_include_files(stdint.h HAVE_STDINT_H)
check_include_files(stdlib.h HAVE_STDLIB_H)
check_function_exists(strchr HAVE_STRCHR)
check_function_exists(strcpy HAVE_STRCPY)
check_include_files(strings.h HAVE_STRINGS_H)
check_include_files(string.h HAVE_STRING_H)
check_function_exists(strstr HAVE_STRSTR)
check_include_files(sys/stat.h HAVE_SYS_STAT_H)
check_include_files(sys/time.h HAVE_SYS_TIME_H)
check_include_files(sys/types.h HAVE_SYS_TYPES_H)
check_include_files(unistd.h HAVE_UNISTD_H)

configure_file(${PROJECT_SOURCE_DIR}/config.h.in.cmake ${PROJECT_SOURCE_DIR}/config.h)
configure_file(${PROJECT_SOURCE_DIR}/polylib-config.in.cmake ${PROJECT_BINARY_DIR}/polylib-config)
configure_file(include/polyVersion.h.in ${PROJECT_BINARY_DIR}/include/polyVersion.h)


# Make_sub_dirs
add_subdirectory(src)
if(USE_NPATCH_LIB)
    add_subdirectory(tools)
endif()
add_subdirectory(doc)

if(enable_example STREQUAL yes)
	add_subdirectory(examples)
endif()

if(enable_test STREQUAL yes)
	add_subdirectory(tests)
endif()


# Install files
install(FILES ${PROJECT_BINARY_DIR}/polylib-config DESTINATION ${PL}/bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )
install(FILES  AUTHORS COPYING ChangeLog README LICENSE DESTINATION ${PL}/share)
