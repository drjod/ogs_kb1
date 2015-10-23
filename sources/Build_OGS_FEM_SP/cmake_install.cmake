# Install script for directory: F:/testingEnvironment/amak/ogs/ogs_kb1/sources

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/OGS")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("F:/testingEnvironment/amak/ogs/ogs_kb1/sources/Build_OGS_FEM_SP/ThirdParty/cmake_install.cmake")
  INCLUDE("F:/testingEnvironment/amak/ogs/ogs_kb1/sources/Build_OGS_FEM_SP/Base/cmake_install.cmake")
  INCLUDE("F:/testingEnvironment/amak/ogs/ogs_kb1/sources/Build_OGS_FEM_SP/MathLib/cmake_install.cmake")
  INCLUDE("F:/testingEnvironment/amak/ogs/ogs_kb1/sources/Build_OGS_FEM_SP/GEO/cmake_install.cmake")
  INCLUDE("F:/testingEnvironment/amak/ogs/ogs_kb1/sources/Build_OGS_FEM_SP/MSH/cmake_install.cmake")
  INCLUDE("F:/testingEnvironment/amak/ogs/ogs_kb1/sources/Build_OGS_FEM_SP/MSHGEOTOOLS/cmake_install.cmake")
  INCLUDE("F:/testingEnvironment/amak/ogs/ogs_kb1/sources/Build_OGS_FEM_SP/FEM/cmake_install.cmake")
  INCLUDE("F:/testingEnvironment/amak/ogs/ogs_kb1/sources/Build_OGS_FEM_SP/GCC/cmake_install.cmake")
  INCLUDE("F:/testingEnvironment/amak/ogs/ogs_kb1/sources/Build_OGS_FEM_SP/FileIO/cmake_install.cmake")
  INCLUDE("F:/testingEnvironment/amak/ogs/ogs_kb1/sources/Build_OGS_FEM_SP/OGS/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "F:/testingEnvironment/amak/ogs/ogs_kb1/sources/Build_OGS_FEM_SP/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "F:/testingEnvironment/amak/ogs/ogs_kb1/sources/Build_OGS_FEM_SP/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
