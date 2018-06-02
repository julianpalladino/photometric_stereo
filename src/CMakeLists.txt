cmake_minimum_required(VERSION 3.5)
project(sstereo)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(3DBUILD_SRC
        3dbuild.cpp
        ppmloader.cpp
        ppmloader.h
        image.cpp
        image.h
        matrix.h
        common.cpp
        common.h
        )

set(CALIBRATE_SRC
        calibrate.cpp
        ppmloader.cpp
        ppmloader.h
        image.cpp
        image.h
        matrix.h
        common.cpp
        common.h
        )

set(TESTS_SRC
        tests.cpp
        matrix.h
        common.cpp
        common.h
        )

# descomentar para debuggear
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")

add_executable(3dbuild ${3DBUILD_SRC})
add_executable(calibrate ${CALIBRATE_SRC})
add_executable(run_tests ${TESTS_SRC})
