cmake_minimum_required(VERSION 3.10.2)
project(MorphometryExamples)

find_package(ITK REQUIRED COMPONENTS Morphometry ITKIOImageBase ITKIONIFTI ITKIONRRD ITKThresholding ITKDistanceMap ITKMesh ITKSmoothing)
include(${ITK_USE_FILE})

add_executable(cleanSurface cleanSurface.cxx)
target_link_libraries(cleanSurface ${ITK_LIBRARIES})

add_executable(highOrderSDT highOrderSDT.cxx)
target_link_libraries(highOrderSDT ${ITK_LIBRARIES})

add_executable(correction correction.cxx)
target_link_libraries(correction ${ITK_LIBRARIES})

add_executable(binaryMorphometry binaryMorphometry.cxx)
target_link_libraries(binaryMorphometry ${ITK_LIBRARIES})

add_executable(quantError quantError.cxx)
target_link_libraries(quantError ${ITK_LIBRARIES})

add_executable(smoothBinary smoothBinary.cxx)
target_link_libraries(smoothBinary ${ITK_LIBRARIES})

add_executable(convergence convergence.cxx)
target_link_libraries(convergence ${ITK_LIBRARIES})

add_executable(noise noise.cxx)
target_link_libraries(noise ${ITK_LIBRARIES})

add_executable(init init.cxx)
target_link_libraries(init ${ITK_LIBRARIES})

add_executable(combine combine.cxx)
target_link_libraries(combine ${ITK_LIBRARIES})

add_executable(narrow_band_order narrow_band_order.cxx)
target_link_libraries(narrow_band_order ${ITK_LIBRARIES})

add_executable(sweeping_order sweeping_order.cxx)
target_link_libraries(sweeping_order ${ITK_LIBRARIES})

add_executable(hsdt hsdt.cxx)
target_link_libraries(hsdt ${ITK_LIBRARIES})

add_executable(rho rho.cxx)
target_link_libraries(rho ${ITK_LIBRARIES})

set(INSTALL_RUNTIME_DESTINATION bin CACHE STRING "Install destination")

install(
  TARGETS cleanSurface highOrderSDT correction binaryMorphometry quantError smoothBinary convergence noise init narrow_band_order sweeping_order hsdt rho
  DESTINATION ${INSTALL_RUNTIME_DESTINATION}
  )
