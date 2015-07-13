::: The Visual Studio project file will be created in
::: /sources/Build/OGS-FEM-5.sln
set USE_VC_EXPRESS=FALSE

::: Create build directory :::
rd /S /Q sources\Build_OGS_FEM_SP
mkdir sources\Build_OGS_FEM_SP

cd sources\Build_OGS_FEM_SP
cmake .. -DOGS_FEM_SP=ON %CMAKE_GENERATOR%
cmake ..

PAUSE