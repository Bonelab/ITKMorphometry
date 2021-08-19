#include <iostream>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkHighOrderSignedDistanceTransformImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

int main(int argc, char * argv[])
{
  if( argc < 4 )
  {
    std::cerr << "Usage: "<< std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFileName> <OutputMeasure> <Threshold>";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  /* Read input Parameters */
  std::string inputFileName = argv[1];
  std::string outputFileName = argv[2];
  int threshold = std::stoi(argv[3]);

  std::cout << "Read in the following parameters:" << std::endl;
  std::cout << "  Input Filename:  " << inputFileName << std::endl;
  std::cout << "  Output Fileanme: " << outputFileName << std::endl;
  std::cout << "  Threshold:       " << threshold << std::endl;
  std::cout << std::endl;

  /* Setup Types */
  constexpr unsigned int ImageDimension = 3;
  using InputPixelType = unsigned int;
  using InputImageType = itk::Image<InputPixelType, ImageDimension>;
  using OutputPixelType = double;
  using OutputImageType = itk::Image<OutputPixelType, ImageDimension>;

  using ReaderType = itk::ImageFileReader< InputImageType >;
  using WriterType = itk::ImageFileWriter< OutputImageType >;
  using BinaryThresholdImageFilterType = itk::BinaryThresholdImageFilter< InputImageType, InputImageType >;
  using SDTType = itk::HighOrderSignedDistanceTransformImageFilter< InputImageType, OutputImageType >;

  /* Read */
  std::cout << "Reading in " << inputFileName << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputFileName);
  reader->Update();

  /* Threshold */
  std::cout << "Thresholding at " << threshold << std::endl;
  BinaryThresholdImageFilterType::Pointer thresh = BinaryThresholdImageFilterType::New();
  thresh->SetInput(reader->GetOutput());
  thresh->SetLowerThreshold(threshold);
  thresh->SetUpperThreshold(threshold);
  thresh->SetInsideValue(1);
  thresh->SetOutsideValue(0);
  thresh->Update();

  /* SDT */
  std::cout << "High order SDT" << std::endl;
  SDTType::Pointer fitler = SDTType::New();
  fitler->SetInput(thresh->GetOutput());
  fitler->Update();

  /* Write */
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(fitler->GetOutput());
  writer->SetFileName(outputFileName);

  std::cout << "Writing results to " << outputFileName << std::endl;
  writer->Write();

  return EXIT_SUCCESS;
}
