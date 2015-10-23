/**
 * \file BuildInfo.h.in
 * 24/08/2010 LB Initial implementation
 * #defines which gets set through CMake
 */

#ifndef BUILDINFO_H
#define BUILDINFO_H

/* #undef CMAKE_CMD_ARGS */
/* #undef SVN_REVISION */
/* #undef GIT_COMMIT_INFO */
/* #undef GIT_BRANCH_INFO */
#define BUILD_TIMESTAMP "2015-10-23"
/* #undef CMAKE_RUN_TIMESTAMP */
#define CMAKE_SYSTEM "Windows-6.1-x32"
#define CMAKE_SYSTEM_PROCESSOR "AMD64"
#define CMAKE_CXX_COMPILER "C:/Program Files (x86)/Microsoft Visual Studio 12.0/VC/bin/cl.exe"
/* #undef GCC_VERSION */
#define CMAKE_GENERATOR "Visual Studio 12"
/* #undef CMAKE_BUILD_TYPE */
#define CMAKE_CXX_FLAGS " /DWIN32 /D_WINDOWS /W3 /GR /EHsc /W3 /wd4290 /wd4267"
#define CMAKE_CXX_FLAGS_RELEASE "/MD /O2 /Ob2 /D NDEBUG /MP"
#define CMAKE_CXX_FLAGS_DEBUG "/D_DEBUG /MDd /Zi /Ob0 /Od   /MP"

#endif // BUILDINFO_H
