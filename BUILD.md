------------------
Build instructions
------------------

Required tools:
* [CMake]
* [git]
* [GCC] //g++ g++4.8

sudo apt-get install cmake
sudo apt-get install git
sudo apt-get install g++

-----------------
Linux compilation
-----------------

[Ubuntu]

```
## make directory
mkdir i23d
cd i23d
`cp i23dMVS i23dSFM ceres-solver PoissonReconTrim plySimp to i23d`
##
#Setup the required external library:
##

#Ceres (Required)
sudo apt-get install libgoogle-glog-dev libatlas-base-dev libeigen3-dev libsuitesparse-dev
#git clone https://ceres-solver.googlesource.com/ceres-solver
mkdir ceres_build
cd ceres_build/
cmake . ../ceres-solver/
make -j
sudo make install
cd ..

#CGAL (Required)
sudo apt-get install libcgal-dev

#OpenCV (Required)
sudo apt-get install libopencv-dev

#VCGLib (Not Required)
#sudo apt-get install subversion
#svn checkout svn://svn.code.sf.net/p/vcg/code/trunk/vcglib vcglib

#MESA (Required)
sudo apt-get install mesa-common-dev

#Boost (Required)
sudo apt-get install libboost-iostreams-dev libboost-program-options-dev libboost-system-dev libboost-serialization-dev

#I23DSFM (Required)
sudo apt-get install libpng-dev libjpeg-dev libtiff-dev libxxf86vm1 libxxf86vm-dev libxi-dev libxrandr-dev graphviz
current_path=`pwd`
mkdir i23dSFM_build
cd i23dSFM_build
cmake -DCMAKE_BUILD_TYPE=RELEASE . ../i23dSFM/ -DCMAKE_INSTALL_PREFIX=$current_path/i23dSFM_install
make -j
make install
cd ..

#I23DMVS build
sudo ln -s /usr/lib/x86_64-linux-gnu/libGLU.so.1 /usr/lib/x86_64-linux-gnu/libGLU.so
sudo ln -s /usr/lib/x86_64-linux-gnu/mesa/libGL.so.1 /usr/lib/x86_64-linux-gnu/libGL.so

mkdir i23dMVS_build
cd i23dMVS_build
`remove include_install in i23dSFM_install/share/i23dSFM/cmake/I23dSFMConfig.cmake`
cmake . ../i23dMVS -DCMAKE_BUILD_TYPE=RELEASE -DVCG_DIR="../i23dMVS/third_party/vcglib" -DOpenCV_CAN_BREAK_BINARY_COMPATIBILITY=OFF -DI23dSFM_DIR:STRING="$current_path/i23dSFM_install/share/i23dSFM/cmake/"
make -j
cd ..

cd PoissonReconTrim
make clean
make
cd ..

