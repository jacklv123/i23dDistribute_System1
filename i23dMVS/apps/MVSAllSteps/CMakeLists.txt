FIND_PACKAGE(I23dSFM QUIET)
IF (I23DSFM_FOUND)
	INCLUDE_DIRECTORIES(${I23DSFM_INCLUDE_DIRS})
	add_definitions(-D_USE_I23DSFM)
	set(LIBS_DEPEND "MVS;${I23DSFM_LIBRARIES}")
ELSE()
	set(LIBS_DEPEND "MVS")
	MESSAGE("I23DSFM_NOT FOUND : I23dSFM importer with JSON support will not be build")
ENDIF()

if(MSVC)
	FILE(GLOB LIBRARY_FILES_C "*.cpp" "*.rc")
else()
	FILE(GLOB LIBRARY_FILES_C "*.cpp")
endif()

FILE(GLOB LIBRARY_FILES_H "*.h" "*.inl")
#SET(LIBRARY_FILES_C MVSAllSteps.cpp)
#SET(LIBRARY_FILES_H MVSAllSteps.h)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#cxx_executable_with_flags_no_pch(MVSAllSteps "Apps" "${cxx_default}" "${LIBS_DEPEND}"  ${LIBRARY_FILES_C} ${LIBRARY_FILES_H})
ADD_LIBRARY(MVSAllSteps ${LIBRARY_FILES_C} ${LIBRARY_FILES_H})
#cxx_library_with_type_no_pch(MVSAllSteps "Apps" "" "${cxx_default}"  ${LIBRARY_FILES_C} ${LIBRARY_FILES_H})

target_link_libraries(MVSAllSteps ${LIBS_DEPEND})


if(I23D_USE_BREAKPAD AND BREAKPAD_FOUND)
	target_link_libraries(MVSAllSteps ${BREAKPAD_LIBS})
endif()
