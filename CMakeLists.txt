cmake_minimum_required(VERSION 3.6)
project(vchsm VERSION 0.1 LANGUAGES CXX)

set(SRC_DIR CppAlgo/src)

# SIMD support (AVX here) (for all build types) and fast floating-point math mode
# this setting applies to all the targets defined after this
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
    add_compile_options(/arch:AVX /fp:fast /Oi /Ot /GL)
else()
    #set(CCMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
    add_compile_options(-mavx -ffast-math -Ofast)
endif()

add_executable(vchsm
	${SRC_DIR}/aaalsf.cpp
	${SRC_DIR}/automaticth2wfw.cpp
	${SRC_DIR}/bylevdurb.cpp
	${SRC_DIR}/bymcaquat.cpp
	${SRC_DIR}/clustkmeans.cpp
	${SRC_DIR}/comp_delta.cpp
	${SRC_DIR}/conv.cpp
	${SRC_DIR}/convert_C.cpp
	${SRC_DIR}/decomposephase.cpp
	${SRC_DIR}/detalsf.cpp
    ${SRC_DIR}/DynamicTimeWarping.cpp
	${SRC_DIR}/extractmfcc.cpp
	${SRC_DIR}/f0analysisbyboersma.cpp
	${SRC_DIR}/harmonicanalysis.cpp
	${SRC_DIR}/HSManalyze.cpp
	${SRC_DIR}/HSManalyze_mfcc.cpp
	${SRC_DIR}/HSMptraining.cpp
	${SRC_DIR}/HSMsynthesize.cpp
	${SRC_DIR}/HSMwfwconvert.cpp
	${SRC_DIR}/linphaseterm.cpp
	${SRC_DIR}/lsfadap.cpp
	CppAlgo/example/main.cpp
	${SRC_DIR}/mfcc.cpp
	${SRC_DIR}/modelSerialization.cpp
	${SRC_DIR}/polarityanalysis.cpp
	${SRC_DIR}/polytseb.cpp
	${SRC_DIR}/PrepareParallelData.cpp
	${SRC_DIR}/RealTimeConverter.cpp
	${SRC_DIR}/stochalsf.cpp
	${SRC_DIR}/stochasticanalysis.cpp
	${SRC_DIR}/sumcos.cpp
	${SRC_DIR}/train_C.cpp
	${SRC_DIR}/trifbank.cpp
	${SRC_DIR}/vec2frames.cpp
	${SRC_DIR}/WavRW.cpp
)

# C++ 11
target_compile_features(vchsm PUBLIC cxx_std_14)
set_target_properties(vchsm PROPERTIES CXX_EXTENSIONS OFF)

# include directories
target_include_directories(vchsm PRIVATE CppAlgo/external/Eigen3.3.4)
target_include_directories(vchsm PRIVATE CppAlgo/include/vchsm)

# macro NDEBUG in release mode
target_compile_definitions(vchsm PRIVATE "$<$<CONFIG:RELEASE>:NDEBUG>")

# link to the thread lib for <thread>
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)
target_link_libraries(vchsm Threads::Threads)



# set the build type to default release if not specified by the user
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
	message(STATUS "Setting build type to ${CMAKE_BUILD_TYPE} since none was specified")
  endif()
endif()
