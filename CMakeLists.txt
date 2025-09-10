cmake_minimum_required(VERSION 3.10)
project(hysteresis_analyzer CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# source files
set(SOURCES
    src/main.cpp
    src/DataLoader.cpp
    src/HysteresisAnalyzer.cpp
    src/ParametricModel.cpp
    src/ModelFitter.cpp
    src/LMOptimizer.cpp
    src/PolynomialModel.cpp
   )

# header files
set(HEADERS
    include/DataLoader.h
    include/HysteresisAnalyzer.h
    include/ParametricModel.h
    include/ModelFitter.h
    include/LMOptimizer.h
    include/PolynomialModel.h
   )

include_directories(include)

add_executable(hysteresis_analyzer ${SOURCES} ${HEADERS})

# link the std++fs excplicitly for GCC < 9
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    target_link_libraries(hysteresis_analyzer stdc++fs)
endif()

# enable testing for later
enable_testing()
