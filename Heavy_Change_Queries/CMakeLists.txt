cmake_minimum_required (VERSION 2.6)
project(StreamClassifier)

if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "GCC detected, adding compile flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra -Wredundant-decls")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2 -mssse3 -msse4.1 -msse4.2 -mavx -mbmi -march=native")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g")
    ADD_DEFINITIONS(-DUNIX)
endif(CMAKE_COMPILER_IS_GNUCXX)

if (UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
endif(UNIX)

ADD_DEFINITIONS(-DPROJECT_ROOT="${CMAKE_SOURCE_DIR}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(SOURCE_FILES src/main.cpp)


add_executable(demo src/main.cpp src/Csketch.h src/Space_allocation.h src/Pbsketch_distr_space.h src/CU_For_SC.h src/FlowRadar_without_packetcount.h src/XY_FlowRadar.h src/Heavy_guardian_flowradar.h src/HG_XY_heavy.h src/HG_for_XY.h src/Heavy_guardian_list.h src/HG_for_XY_without_auxiliary.h "src/HG(heavy_light).h")

