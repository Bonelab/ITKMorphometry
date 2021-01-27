cmake_minimum_required(VERSION 3.10.2)
project(MorphometryExamples)

find_package(ITK REQUIRED COMPONENTS Morphometry ITKIOImageBase ITKIONIFTI ITKIONRRD ITKThresholding)
include(${ITK_USE_FILE})

add_executable(cleanSurface cleanSurface.cxx)
target_link_libraries(cleanSurface ${ITK_LIBRARIES})

set(INSTALL_RUNTIME_DESTINATION bin CACHE STRING "Install destination")

install(
  TARGETS cleanSurface
  DESTINATION ${INSTALL_RUNTIME_DESTINATION}
  )