::: The Visual Studio project file will be created in
::: /Build/OGS-FEM-5.sln
set USE_VC_EXPRESS=FALSE

::: Create build directory :::
rd /S /Q Build_OGS_FEM
mkdir Build_OGS_FEM

cd Build_OGS_FEM
cmake ../sources -DOGS_FEM=ON %CMAKE_GENERATOR%

PAUSE
