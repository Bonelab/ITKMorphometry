cmake_minimum_required(VERSION 3.10.2)
project(Morphometry)

set(Morphometry_LIBRARIES Morphometry)

if(NOT ITK_SOURCE_DIR)
  find_package(ITK REQUIRED ITKThresholding ITKDistanceMap)
  list(APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR})
  include(ITKModuleExternal)
else()
  set(ITK_DIR ${CMAKE_BINARY_DIR})
  itk_module_impl()
endif()
