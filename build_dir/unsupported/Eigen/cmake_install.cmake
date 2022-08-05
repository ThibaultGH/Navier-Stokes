# Install script for directory: /home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/AdolcForward"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/AlignedVector3"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/ArpackSupport"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/AutoDiff"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/BVH"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/EulerAngles"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/FFT"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/IterativeSolvers"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/KroneckerProduct"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/LevenbergMarquardt"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/MatrixFunctions"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/MoreVectorization"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/MPRealSupport"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/NonLinearOptimization"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/NumericalDiff"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/OpenGLSupport"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/Polynomials"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/Skyline"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/SparseExtra"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/SpecialFunctions"
    "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/Splines"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/source_dir/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/thibault/Dropbox/M2 - Mathématiques/Cours Fondamentaux/Resolution EDP par EF/TP4/build_dir/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

