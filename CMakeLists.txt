cmake_minimum_required (VERSION 3.8)
project (rt)

# get rid of annoying MSVC warnings.
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

set (CMAKE_CXX_STANDARD 11)

add_executable(rt main.cpp Bitmap.cpp my_math.h Primitives.cpp Primitives.h Model.cpp Model.h)

target_link_libraries(rt ${ALL_LIBS} )

set (CMAKE_CXX_FLAGS "-fopenmp -O3")
