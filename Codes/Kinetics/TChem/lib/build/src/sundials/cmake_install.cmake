# Install script for directory: /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/src/sundials

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0-install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
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

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  MESSAGE("
Install shared components
")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sundials" TYPE FILE FILES
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/sundials/sundials_band.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/sundials/sundials_dense.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/sundials/sundials_direct.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/sundials/sundials_iterative.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/sundials/sundials_math.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/sundials/sundials_nvector.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/sundials/sundials_fnvector.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/sundials/sundials_spbcgs.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/sundials/sundials_spgmr.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/sundials/sundials_sptfqmr.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/sundials/sundials_types.h"
    )
endif()

