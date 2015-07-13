cd $SOURCE_LOCATION/..
rm -rf build_fem
rm -rf build_gui
mkdir -vp build_fem
mkdir -vp build_gui

cd build_fem
cmake -DOGS_FEM=ON -G "$CMAKE_GENERATOR" $SOURCE_LOCATION
cmake $SOURCE_LOCATION
#if [ "$OSTYPE" == 'msys' ]; then
#	if [ -f OGS.sln ]; then
#		devenv OGS.sln &
#	else
#		echo "CMake configuration of FEM project went wrong. Aborting..."
#		exit 1
#	fi
#fi

cd $SOURCE_LOCATION/../

cd build_gui
cmake -DOGS_USE_QT=ON -G "$CMAKE_GENERATOR" $SOURCE_LOCATION
cmake $SOURCE_LOCATION
if [ "$OSTYPE" == 'msys' ]; then
	if [ -f OGS.sln ]; then
		devenv OGS.sln &
	else
		echo "CMake configuration of GUI project went wrong. Aborting..."
		exit 1
	fi
fi
