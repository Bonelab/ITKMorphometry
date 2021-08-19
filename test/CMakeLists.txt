itk_module_test()

set(MorphometryTests
  itkMinimalStandardRandomVariateGeneratorTest.cxx
  itkMeanCurvatureImageFilterTest.cxx
  )

CreateTestDriver(Morphometry "${Morphometry-Test_LIBRARIES}" "${MorphometryTests}")

itk_add_test(NAME itkMinimalStandardRandomVariateGeneratorTest
  COMMAND MorphometryTestDriver itkMinimalStandardRandomVariateGeneratorTest
  )

itk_add_test(NAME itkMeanCurvatureImageFilterTest
  COMMAND MorphometryTestDriver
    --compare
    DATA{Baseline/itkMeanCurvatureImageFilterTestOutput.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkMeanCurvatureImageFilterTestOutput.mha
  itkMeanCurvatureImageFilterTest
    ${ITK_TEST_OUTPUT_DIR}/itkMeanCurvatureImageFilterTestOutput.mha
  )

set(MorphometryUnitTests
    itkFiniteDifferenceUnitTest.cxx
    itkSignalFunctionsUnitTest.cxx
    itkHeavisideImageFilterUnitTest.cxx
    itkDiracDeltaImageFilterUnitTest.cxx
    itkCleanBinarySurfaceImageFilterUnitTest.cxx
  )

CreateGoogleTestDriver(MorphometryUnitTests "${Morphometry-Test_LIBRARIES}" "${MorphometryUnitTests}")
