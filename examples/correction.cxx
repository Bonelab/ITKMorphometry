#include <iostream>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkHighOrderSignedDistanceTransformCorrectionImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkInitializeSignedDistanceTransform.h"
#include "itkHighOrderSignedFastSweeping.hxx"

int main(int argc, char * argv[])
{
  if( argc < 5 )
  {
    std::cerr << "Usage: "<< std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFileName> <OutputMeasure> <Threshold> <NarrowbandSize>";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  /* Read input Parameters */
  std::string inputFileName = argv[1];
  std::string outputFileName = argv[2];
  int threshold = std::stoi(argv[3]);
  double narrowbandSize = std::stof(argv[4]);

  std::cout << "Read in the following parameters:" << std::endl;
  std::cout << "  Input Filename:  " << inputFileName << std::endl;
  std::cout << "  Output Fileanme: " << outputFileName << std::endl;
  std::cout << "  Threshold:       " << threshold << std::endl;
  std::cout << "  Narrowband size: " << narrowbandSize << std::endl;
  std::cout << std::endl;

  /* Setup Types */
  constexpr unsigned int ImageDimension = 3;
  using InputPixelType = unsigned int;
  using InputImageType = itk::Image<InputPixelType, ImageDimension>;
  using OutputPixelType = double;
  using OutputImageType = itk::Image<OutputPixelType, ImageDimension>;

  using ReaderType = itk::ImageFileReader< InputImageType >;
  using WriterType = itk::ImageFileWriter< OutputImageType >;
  using SignedMaurerType = itk::InitializeSignedDistanceTransform< InputImageType, OutputImageType >;
  // using SDTType = itk::HighOrderSignedDistanceTransformCorrectionImageFilter< OutputImageType >;
  using SDTType = itk::HighOrderSignedFastSweeping< OutputImageType >;

  /* Read */
  std::cout << "Reading in " << inputFileName << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputFileName);
  reader->Update();

  /* SDT */
  std::cout << "Signed Distance Transform" << std::endl;
  SignedMaurerType::Pointer sdt = SignedMaurerType::New();
  sdt->SetInput(reader->GetOutput());
  sdt->SetLabel(threshold);
  sdt->SetDither(false);
  sdt->Update();

  /* Correction */
  std::cout << "Correct SDT" << std::endl;
  SDTType::Pointer fitler = SDTType::New();
  fitler->SetInput(sdt->GetOutput());
  fitler->SetNarrowbandSize(narrowbandSize);
  fitler->Update();

  /* Write */
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(fitler->GetOutput());
  // writer->SetInput(sdt->GetOutput());
  writer->SetFileName(outputFileName);

  std::cout << "Writing results to " << outputFileName << std::endl;
  writer->Write();

  return EXIT_SUCCESS;
}
