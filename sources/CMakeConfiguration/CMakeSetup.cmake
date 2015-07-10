
# Set additional CMake modules path
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeConfiguration")
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Load addional modules
INCLUDE(UseBackportedModules)
INCLUDE(OptionRequires)
INCLUDE(CppcheckTargets)

# Adds useful macros and variables
INCLUDE( CMakeConfiguration/Macros.cmake )


# Suppress warning on setting policies
CMAKE_POLICY(SET CMP0011 OLD)

# Suppress warning on add_subdirectory(dir) where dir contains no CMakeLists.txt
IF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 2.7)
	CMAKE_POLICY(SET CMP0014 OLD)
ENDIF ()

# Provide a way for Visual Studio Express users to turn OFF the new FOLDER
# organization feature. Default to ON for non-Express users. Express users must
# explicitly turn off this option to build CMake in the Express IDE...
OPTION(CMAKE_USE_FOLDERS "Enable folder grouping of projects in IDEs." ON)
MARK_AS_ADVANCED(CMAKE_USE_FOLDERS)
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ${CMAKE_USE_FOLDERS})