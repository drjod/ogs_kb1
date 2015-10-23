/**
 * \file Configure.h.in
 *
 * #defines which gets set through CMake
 */
 #ifndef CONFIGURE_H
 #define CONFIGURE_H

#define OGS_QT_VERSION 
#define SOURCEPATH "F:/testingEnvironment/amak/ogs/ogs_kb1/sources"
#define BUILDPATH "F:/testingEnvironment/amak/ogs/ogs_kb1/sources/Build_OGS_FEM"
#define TESTDATAPATH "TESTDATA_DIR_FOUND-NOTFOUND"

#define OGS_VERSION "5.6(CL/TN)"
#define OGS_DATE "07.07.2015"
/* #undef QT_USE_QTXMLPATTERNS */

// for tests
#define OGS_EXECUTABLE "F:/testingEnvironment/amak/ogs/ogs_kb1/sources/Build_OGS_FEM/bin/release/ogs"
#define PUT_TMP_DIR_IN "F:/testingEnvironment/amak/ogs/ogs_kb1/sources/Build_OGS_FEM/tests/"
/* #undef JENKINS_URL */
/* #undef JENKINS_JOB_NAME */
#define PROCESSOR_COUNT 12

#endif // CONFIGURE_H
