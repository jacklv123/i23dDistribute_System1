# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/gaomz/CudaSift-Maxwell

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/gaomz/CudaSift-Maxwell/build

# Include any dependencies generated for this target.
include CMakeFiles/extractsift.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/extractsift.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/extractsift.dir/flags.make

CMakeFiles/cudasift.dir/./cudasift_generated_cudaImage.cu.o: CMakeFiles/cudasift.dir/cudasift_generated_cudaImage.cu.o.depend
CMakeFiles/cudasift.dir/./cudasift_generated_cudaImage.cu.o: CMakeFiles/cudasift.dir/cudasift_generated_cudaImage.cu.o.cmake
CMakeFiles/cudasift.dir/./cudasift_generated_cudaImage.cu.o: ../cudaImage.cu
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gaomz/CudaSift-Maxwell/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Building NVCC (Device) object CMakeFiles/cudasift.dir//./cudasift_generated_cudaImage.cu.o"
	cd /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir && /usr/bin/cmake -E make_directory /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir//.
	cd /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir && /usr/bin/cmake -D verbose:BOOL=$(VERBOSE) -D build_configuration:STRING= -D generated_file:STRING=/home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir//./cudasift_generated_cudaImage.cu.o -D generated_cubin_file:STRING=/home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir//./cudasift_generated_cudaImage.cu.o.cubin.txt -P /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir//cudasift_generated_cudaImage.cu.o.cmake

CMakeFiles/cudasift.dir/./cudasift_generated_cudaSiftH.cu.o: CMakeFiles/cudasift.dir/cudasift_generated_cudaSiftH.cu.o.depend
CMakeFiles/cudasift.dir/./cudasift_generated_cudaSiftH.cu.o: CMakeFiles/cudasift.dir/cudasift_generated_cudaSiftH.cu.o.cmake
CMakeFiles/cudasift.dir/./cudasift_generated_cudaSiftH.cu.o: ../cudaSiftH.cu
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gaomz/CudaSift-Maxwell/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Building NVCC (Device) object CMakeFiles/cudasift.dir//./cudasift_generated_cudaSiftH.cu.o"
	cd /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir && /usr/bin/cmake -E make_directory /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir//.
	cd /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir && /usr/bin/cmake -D verbose:BOOL=$(VERBOSE) -D build_configuration:STRING= -D generated_file:STRING=/home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir//./cudasift_generated_cudaSiftH.cu.o -D generated_cubin_file:STRING=/home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir//./cudasift_generated_cudaSiftH.cu.o.cubin.txt -P /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir//cudasift_generated_cudaSiftH.cu.o.cmake

CMakeFiles/cudasift.dir/./cudasift_generated_matching.cu.o: CMakeFiles/cudasift.dir/cudasift_generated_matching.cu.o.depend
CMakeFiles/cudasift.dir/./cudasift_generated_matching.cu.o: CMakeFiles/cudasift.dir/cudasift_generated_matching.cu.o.cmake
CMakeFiles/cudasift.dir/./cudasift_generated_matching.cu.o: ../matching.cu
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gaomz/CudaSift-Maxwell/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Building NVCC (Device) object CMakeFiles/cudasift.dir//./cudasift_generated_matching.cu.o"
	cd /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir && /usr/bin/cmake -E make_directory /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir//.
	cd /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir && /usr/bin/cmake -D verbose:BOOL=$(VERBOSE) -D build_configuration:STRING= -D generated_file:STRING=/home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir//./cudasift_generated_matching.cu.o -D generated_cubin_file:STRING=/home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir//./cudasift_generated_matching.cu.o.cubin.txt -P /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/cudasift.dir//cudasift_generated_matching.cu.o.cmake

CMakeFiles/extractsift.dir/./extractsift_generated_cudaImage.cu.o: CMakeFiles/extractsift.dir/extractsift_generated_cudaImage.cu.o.depend
CMakeFiles/extractsift.dir/./extractsift_generated_cudaImage.cu.o: CMakeFiles/extractsift.dir/extractsift_generated_cudaImage.cu.o.cmake
CMakeFiles/extractsift.dir/./extractsift_generated_cudaImage.cu.o: ../cudaImage.cu
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gaomz/CudaSift-Maxwell/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Building NVCC (Device) object CMakeFiles/extractsift.dir//./extractsift_generated_cudaImage.cu.o"
	cd /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir && /usr/bin/cmake -E make_directory /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir//.
	cd /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir && /usr/bin/cmake -D verbose:BOOL=$(VERBOSE) -D build_configuration:STRING= -D generated_file:STRING=/home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir//./extractsift_generated_cudaImage.cu.o -D generated_cubin_file:STRING=/home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir//./extractsift_generated_cudaImage.cu.o.cubin.txt -P /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir//extractsift_generated_cudaImage.cu.o.cmake

CMakeFiles/extractsift.dir/./extractsift_generated_cudaSiftH.cu.o: CMakeFiles/extractsift.dir/extractsift_generated_cudaSiftH.cu.o.depend
CMakeFiles/extractsift.dir/./extractsift_generated_cudaSiftH.cu.o: CMakeFiles/extractsift.dir/extractsift_generated_cudaSiftH.cu.o.cmake
CMakeFiles/extractsift.dir/./extractsift_generated_cudaSiftH.cu.o: ../cudaSiftH.cu
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gaomz/CudaSift-Maxwell/build/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Building NVCC (Device) object CMakeFiles/extractsift.dir//./extractsift_generated_cudaSiftH.cu.o"
	cd /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir && /usr/bin/cmake -E make_directory /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir//.
	cd /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir && /usr/bin/cmake -D verbose:BOOL=$(VERBOSE) -D build_configuration:STRING= -D generated_file:STRING=/home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir//./extractsift_generated_cudaSiftH.cu.o -D generated_cubin_file:STRING=/home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir//./extractsift_generated_cudaSiftH.cu.o.cubin.txt -P /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir//extractsift_generated_cudaSiftH.cu.o.cmake

CMakeFiles/extractsift.dir/./extractsift_generated_matching.cu.o: CMakeFiles/extractsift.dir/extractsift_generated_matching.cu.o.depend
CMakeFiles/extractsift.dir/./extractsift_generated_matching.cu.o: CMakeFiles/extractsift.dir/extractsift_generated_matching.cu.o.cmake
CMakeFiles/extractsift.dir/./extractsift_generated_matching.cu.o: ../matching.cu
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gaomz/CudaSift-Maxwell/build/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Building NVCC (Device) object CMakeFiles/extractsift.dir//./extractsift_generated_matching.cu.o"
	cd /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir && /usr/bin/cmake -E make_directory /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir//.
	cd /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir && /usr/bin/cmake -D verbose:BOOL=$(VERBOSE) -D build_configuration:STRING= -D generated_file:STRING=/home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir//./extractsift_generated_matching.cu.o -D generated_cubin_file:STRING=/home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir//./extractsift_generated_matching.cu.o.cubin.txt -P /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir//extractsift_generated_matching.cu.o.cmake

CMakeFiles/extractsift.dir/./extractsift_intermediate_link.o: CMakeFiles/extractsift.dir/./extractsift_generated_cudaImage.cu.o
CMakeFiles/extractsift.dir/./extractsift_intermediate_link.o: CMakeFiles/extractsift.dir/./extractsift_generated_cudaSiftH.cu.o
CMakeFiles/extractsift.dir/./extractsift_intermediate_link.o: CMakeFiles/extractsift.dir/./extractsift_generated_matching.cu.o
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gaomz/CudaSift-Maxwell/build/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Building NVCC intermediate link file CMakeFiles/extractsift.dir/./extractsift_intermediate_link.o"
	/usr/local/cuda/bin/nvcc -arch=sm_35 -m64 -ccbin "/usr/bin/cc" -dlink /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir//./extractsift_generated_cudaImage.cu.o /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir//./extractsift_generated_cudaSiftH.cu.o /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir//./extractsift_generated_matching.cu.o -o /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir/./extractsift_intermediate_link.o

CMakeFiles/extractsift.dir/geomFuncs.cpp.o: CMakeFiles/extractsift.dir/flags.make
CMakeFiles/extractsift.dir/geomFuncs.cpp.o: ../geomFuncs.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gaomz/CudaSift-Maxwell/build/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/extractsift.dir/geomFuncs.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/extractsift.dir/geomFuncs.cpp.o -c /home/gaomz/CudaSift-Maxwell/geomFuncs.cpp

CMakeFiles/extractsift.dir/geomFuncs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/extractsift.dir/geomFuncs.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/gaomz/CudaSift-Maxwell/geomFuncs.cpp > CMakeFiles/extractsift.dir/geomFuncs.cpp.i

CMakeFiles/extractsift.dir/geomFuncs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/extractsift.dir/geomFuncs.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/gaomz/CudaSift-Maxwell/geomFuncs.cpp -o CMakeFiles/extractsift.dir/geomFuncs.cpp.s

CMakeFiles/extractsift.dir/geomFuncs.cpp.o.requires:
.PHONY : CMakeFiles/extractsift.dir/geomFuncs.cpp.o.requires

CMakeFiles/extractsift.dir/geomFuncs.cpp.o.provides: CMakeFiles/extractsift.dir/geomFuncs.cpp.o.requires
	$(MAKE) -f CMakeFiles/extractsift.dir/build.make CMakeFiles/extractsift.dir/geomFuncs.cpp.o.provides.build
.PHONY : CMakeFiles/extractsift.dir/geomFuncs.cpp.o.provides

CMakeFiles/extractsift.dir/geomFuncs.cpp.o.provides.build: CMakeFiles/extractsift.dir/geomFuncs.cpp.o

CMakeFiles/extractsift.dir/extractSift.cpp.o: CMakeFiles/extractsift.dir/flags.make
CMakeFiles/extractsift.dir/extractSift.cpp.o: ../extractSift.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gaomz/CudaSift-Maxwell/build/CMakeFiles $(CMAKE_PROGRESS_9)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/extractsift.dir/extractSift.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/extractsift.dir/extractSift.cpp.o -c /home/gaomz/CudaSift-Maxwell/extractSift.cpp

CMakeFiles/extractsift.dir/extractSift.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/extractsift.dir/extractSift.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/gaomz/CudaSift-Maxwell/extractSift.cpp > CMakeFiles/extractsift.dir/extractSift.cpp.i

CMakeFiles/extractsift.dir/extractSift.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/extractsift.dir/extractSift.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/gaomz/CudaSift-Maxwell/extractSift.cpp -o CMakeFiles/extractsift.dir/extractSift.cpp.s

CMakeFiles/extractsift.dir/extractSift.cpp.o.requires:
.PHONY : CMakeFiles/extractsift.dir/extractSift.cpp.o.requires

CMakeFiles/extractsift.dir/extractSift.cpp.o.provides: CMakeFiles/extractsift.dir/extractSift.cpp.o.requires
	$(MAKE) -f CMakeFiles/extractsift.dir/build.make CMakeFiles/extractsift.dir/extractSift.cpp.o.provides.build
.PHONY : CMakeFiles/extractsift.dir/extractSift.cpp.o.provides

CMakeFiles/extractsift.dir/extractSift.cpp.o.provides.build: CMakeFiles/extractsift.dir/extractSift.cpp.o

# Object files for target extractsift
extractsift_OBJECTS = \
"CMakeFiles/extractsift.dir/geomFuncs.cpp.o" \
"CMakeFiles/extractsift.dir/extractSift.cpp.o"

# External object files for target extractsift
extractsift_EXTERNAL_OBJECTS = \
"/home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir/./extractsift_generated_cudaImage.cu.o" \
"/home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir/./extractsift_generated_cudaSiftH.cu.o" \
"/home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir/./extractsift_generated_matching.cu.o" \
"/home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir/./extractsift_intermediate_link.o"

extractsift: CMakeFiles/extractsift.dir/geomFuncs.cpp.o
extractsift: CMakeFiles/extractsift.dir/extractSift.cpp.o
extractsift: CMakeFiles/extractsift.dir/./extractsift_generated_cudaImage.cu.o
extractsift: CMakeFiles/extractsift.dir/./extractsift_generated_cudaSiftH.cu.o
extractsift: CMakeFiles/extractsift.dir/./extractsift_generated_matching.cu.o
extractsift: CMakeFiles/extractsift.dir/./extractsift_intermediate_link.o
extractsift: CMakeFiles/extractsift.dir/build.make
extractsift: /usr/local/cuda/lib64/libcudart.so
extractsift: /usr/local/cuda/lib64/libcudadevrt.a
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_videostab.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_video.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_ts.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_superres.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_stitching.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_photo.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_ocl.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_objdetect.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_ml.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_legacy.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_imgproc.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_highgui.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_gpu.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_flann.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_features2d.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_core.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_contrib.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_calib3d.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_photo.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_legacy.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_video.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_objdetect.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_ml.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_calib3d.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_features2d.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_highgui.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_imgproc.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_flann.so.2.4.8
extractsift: /usr/lib/x86_64-linux-gnu/libopencv_core.so.2.4.8
extractsift: CMakeFiles/extractsift.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable extractsift"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/extractsift.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/extractsift.dir/build: extractsift
.PHONY : CMakeFiles/extractsift.dir/build

CMakeFiles/extractsift.dir/requires: CMakeFiles/extractsift.dir/geomFuncs.cpp.o.requires
CMakeFiles/extractsift.dir/requires: CMakeFiles/extractsift.dir/extractSift.cpp.o.requires
.PHONY : CMakeFiles/extractsift.dir/requires

CMakeFiles/extractsift.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/extractsift.dir/cmake_clean.cmake
.PHONY : CMakeFiles/extractsift.dir/clean

CMakeFiles/extractsift.dir/depend: CMakeFiles/cudasift.dir/./cudasift_generated_cudaImage.cu.o
CMakeFiles/extractsift.dir/depend: CMakeFiles/cudasift.dir/./cudasift_generated_cudaSiftH.cu.o
CMakeFiles/extractsift.dir/depend: CMakeFiles/cudasift.dir/./cudasift_generated_matching.cu.o
CMakeFiles/extractsift.dir/depend: CMakeFiles/extractsift.dir/./extractsift_generated_cudaImage.cu.o
CMakeFiles/extractsift.dir/depend: CMakeFiles/extractsift.dir/./extractsift_generated_cudaSiftH.cu.o
CMakeFiles/extractsift.dir/depend: CMakeFiles/extractsift.dir/./extractsift_generated_matching.cu.o
CMakeFiles/extractsift.dir/depend: CMakeFiles/extractsift.dir/./extractsift_intermediate_link.o
	cd /home/gaomz/CudaSift-Maxwell/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/gaomz/CudaSift-Maxwell /home/gaomz/CudaSift-Maxwell /home/gaomz/CudaSift-Maxwell/build /home/gaomz/CudaSift-Maxwell/build /home/gaomz/CudaSift-Maxwell/build/CMakeFiles/extractsift.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/extractsift.dir/depend

