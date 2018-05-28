-----------------
Linux compilation
-----------------

Setup the required external library.
* sudo apt-get install libpng-dev libjpeg-dev libtiff-dev libxxf86vm1 libxxf86vm-dev libxi-dev libxrandr-dev
If you want see the view graph svg logs
* sudo apt-get install graphviz

 $ tar -xvf i23dSFM.tar.gz
 $ mkdir build
 $ cd build
 $ cmake -DCMAKE_BUILD_TYPE=RELEASE -DI23dSFM_BUILD_TESTS=ON -DI23dSFM_BUILD_EXAMPLES=ON . ../src/

Compile the project
 $ make

For a multi-core compilation (Replace NBcore with the number of threads)
 $ make -j NBcore

For test if build successfully
 $ make test

Usage:
 python i23dSFM_Build/software/Sfm/Sfm_SequentialPipeline.py ${image_dir} ${output_dir}
