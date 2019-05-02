cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(project2B)
SET(VTK_DIR /Users/priyakudva/VTK/VTK-6.3.0)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(project2B MACOSX_BUNDLE project2B)

if(VTK_LIBRARIES)
  target_link_libraries(project2B ${VTK_LIBRARIES})
else()
  target_link_libraries(project2B vtkHybrid)
endif()

