::: The Visual Studio project file will be created in
::: /sources/Build/OGS-FEM-5.sln
set USE_VC_EXPRESS=FALSE

::: Create build directory :::
rd /S /Q sources\Build_OGS_FEM
mkdir sources\Build_OGS_FEM

cd sources\Build_OGS_FEM
cmake .. -DOGS_FEM=ON %CMAKE_GENERATOR%
cmake ..

PAUSE