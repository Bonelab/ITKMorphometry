#include <iostream>
#include <vector>
#include <string>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkInitializeTwoPhaseNarrowBandImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkSelectNarrowbandImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"

constexpr unsigned int ImageDimension = 3;
using PixelType = double;
using ImageType = itk::Image<PixelType, ImageDimension>;
using BinaryPixelType = unsigned int;
using BinaryImageType = itk::Image<BinaryPixelType, ImageDimension>;

using ReaderType = itk::ImageFileReader<ImageType>;
using WriterType = itk::ImageFileWriter<ImageType>;
using BinaryWriterType = itk::ImageFileWriter<BinaryImageType>;
using NBType = itk::SelectNarrowbandImageFilter<ImageType, BinaryImageType>;
using InitType = itk::InitializeTwoPhaseNarrowBandImageFilter<ImageType, BinaryImageType>;
using SubtractType = itk::SubtractImageFilter< ImageType, ImageType, ImageType >;
using SmoothingType = itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType>;

int main(int argc, char * argv[])
{
  if( argc < 5 )
  {
    std::cerr << "Usage: "<< std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFileName> <OutputInitFilename> <OutputMaskFilename> <Threshold>";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  /* Read input Parameters */
  std::string inputFileName = argv[1];
  std::string outputInitFilename = argv[2];
  std::string outputMaskFilename = argv[3];
  double threshold = std::stof(argv[4]);
  double sigma = 0.02;

  std::cout << "Read in the following parameters:" << std::endl;
  std::cout << "  Input Filename:              " << inputFileName << std::endl;
  std::cout << "  Output Initialized Fileanme: " << outputInitFilename << std::endl;
  std::cout << "  Output Mask Fileanme:        " << outputMaskFilename << std::endl;
  std::cout << "  Threshold:                   " << threshold << std::endl;
  std::cout << std::endl;

  std::cout << "Reading in " << inputFileName << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputFileName);
  reader->Update();

  std::cout << "Smoothing with size " << sigma << std::endl;
  SmoothingType::Pointer smoother = SmoothingType::New();
  smoother->SetInput(reader->GetOutput());
  smoother->SetSigma(sigma);
  smoother->Update();

  std::cout << "Shifting by " << threshold << std::endl;
  SubtractType::Pointer sub = SubtractType::New();
  sub->SetConstant1(threshold);
  sub->SetInput2(smoother->GetOutput());
  sub->Update();

  std::cout << "Writing out initalization " << "psi.nii" << std::endl;
  WriterType::Pointer writer1 = WriterType::New();
  writer1->SetInput(sub->GetOutput());
  writer1->SetFileName("psi.nii");
  writer1->Update();

  std::cout << "Selecting narrowband voxels" << std::endl;
  NBType::Pointer nbSelect = NBType::New();
  nbSelect->SetInput(sub->GetOutput());
  nbSelect->SetRadius(3);
  nbSelect->Update();

  std::cout << "Running initializer" << std::endl;
  InitType::Pointer init = InitType::New();
  init->SetGrayInput(sub->GetOutput());
  init->SetNarrowbandInput(nbSelect->GetOutput());
  init->SetMaxIterations(100);
  init->SetConvergenceOrder(3);
  init->Update();

  std::cout << "Writing out initalization " << outputInitFilename << std::endl;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(init->GetInitializedNarrowbandOutput());
  writer->SetFileName(outputInitFilename);
  writer->Update();

  std::cout << "Writing out mask " << outputMaskFilename << std::endl;
  BinaryWriterType::Pointer binaryWriter = BinaryWriterType::New();
  binaryWriter->SetInput(init->GetNarrowbandMaskOutput());
  binaryWriter->SetFileName(outputMaskFilename);
  binaryWriter->Update();

  return EXIT_SUCCESS;
}
