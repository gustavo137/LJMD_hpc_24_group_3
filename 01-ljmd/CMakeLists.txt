cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(LJMD VERSION 1.0)

# Use configure_file to set LJMD version
configure_file(LJMDConfig.h.in LJMDConfig.h @ONLY)

# add the executable
add_executable(ljmd-serial.x src/ljmd.c)
target_compile_options(ljmd-serial.x PUBLIC -Wall -g -O3 -ffast-math -fomit-frame-pointer)
# link math.h
target_link_libraries(ljmd-serial.x m)

# Use target_include_directories to include ${PROJECT_BINARY_DIR}
target_include_directories(ljmd-serial.x PUBLIC "${PROJECT_BINARY_DIR}")
