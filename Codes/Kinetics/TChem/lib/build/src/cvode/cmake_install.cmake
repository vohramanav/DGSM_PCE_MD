# Install script for directory: /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/src/cvode

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
Install CVODE
")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build/src/cvode/libsundials_cvode.a")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/cvode" TYPE FILE FILES
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/cvode/cvode_band.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/cvode/cvode_bandpre.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/cvode/cvode_bbdpre.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/cvode/cvode_dense.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/cvode/cvode_diag.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/cvode/cvode_direct.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/cvode/cvode.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/cvode/cvode_spbcgs.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/cvode/cvode_spgmr.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/cvode/cvode_spils.h"
    "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/include/cvode/cvode_sptfqmr.h"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/cvode" TYPE FILE FILES "/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/src/cvode/cvode_impl.h")
endif()

