itk_module_test()

set(MorphometryUnitTests
    itkCleanBinarySurfaceImageFilterUnitTest.cxx
  )

CreateGoogleTestDriver(MorphometryUnitTests "${Morphometry-Test_LIBRARIES}" "${MorphometryUnitTests}")
