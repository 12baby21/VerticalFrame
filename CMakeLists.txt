cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(Paillier_DJN LANGUAGES CXX)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
# Disable C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
# Require full C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# main.cpp
aux_source_directory(./src SOURCE_DIR)
include_directories(./include)
add_executable(vfl lendingclub.cpp ${SOURCE_DIR})


# link libtorch .a .so
target_link_libraries(vfl gmp)
