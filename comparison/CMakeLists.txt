cmake_minimum_required(VERSION 3.10.2)
project(MorphometryExamples)

find_package(ITK REQUIRED COMPONENTS Morphometry ITKIOImageBase ITKIONIFTI ITKIONRRD ITKThresholding ITKDistanceMap ITKMesh ITKBinaryMathematicalMorphology)
include(${ITK_USE_FILE})

add_executable(binaryCloseWithTime binaryCloseWithTime.cxx)
target_link_libraries(binaryCloseWithTime ${ITK_LIBRARIES})

add_executable(cleanBinaryWithTime cleanBinaryWithTime.cxx)
target_link_libraries(cleanBinaryWithTime ${ITK_LIBRARIES})

add_executable(morphometry morphometry.cxx)
target_link_libraries(morphometry ${ITK_LIBRARIES})

set(INSTALL_RUNTIME_DESTINATION bin CACHE STRING "Install destination")

install(
  TARGETS binaryCloseWithTime cleanBinaryWithTime morphometry
  DESTINATION ${INSTALL_RUNTIME_DESTINATION}
  )