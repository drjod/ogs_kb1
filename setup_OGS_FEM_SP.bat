::: The Visual Studio project file will be created in
::: /Build/OGS-FEM_SP

set USE_VC_EXPRESS=FALSE

::: Create build directory :::
rd /S /Q Build_OGS_FEM_SP
mkdir Build_OGS_FEM_SP

cd Build_OGS_FEM_SP
cmake ../sources -DOGS_FEM_SP=ON %CMAKE_GENERATOR%

PAUSE