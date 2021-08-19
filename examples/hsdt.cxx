#include <iostream>
#include <vector>
#include <string>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkInitializeTwoPhaseNarrowBandImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkSelectNarrowbandImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkHighOrderSignedFastSweeping2.h"

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
using SweepingType = itk::HighOrderSignedFastSweeping2<ImageType, BinaryImageType>;


int main(int argc, char * argv[])
{
  if( argc < 6 )
  {
    std::cerr << "Usage: "<< std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFileName> <OutputFileName> <PsiFileName> <Threshold> <StandardDeviation> <DenseFileName>";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  /* Read input Parameters */
  std::string inputFileName = argv[1];
  std::string outputFileName = argv[2];
  std::string psiFileName = argv[3];
  double threshold = std::stof(argv[4]);
  double sigma = std::stof(argv[5]);
  std::string denseFileName = argv[6];

  std::cout << "Read in the following parameters:" << std::endl;
  std::cout << "  Input Filename:      " << inputFileName << std::endl;
  std::cout << "  Output  Fileanme:    " << outputFileName << std::endl;
  std::cout << "  Psi  Fileanme:       " << psiFileName << std::endl;
  std::cout << "  Threshold:           " << threshold << std::endl;
  std::cout << "  Standard Deviation:  " << sigma << std::endl;
  std::cout << "  Density filename:    " << denseFileName << std::endl;
  std::cout << std::endl;

  std::cout << "Reading in " << inputFileName << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputFileName);
  reader->Update();

  SmoothingType::Pointer smoother = SmoothingType::New();
  std::cout << "Smoothing with size " << sigma << std::endl;
  smoother->SetInput(reader->GetOutput());
  smoother->SetSigma(sigma);
  smoother->Update();

  std::cout << "Writing out dense to " << denseFileName << std::endl;
  WriterType::Pointer writeDense = WriterType::New();
  writeDense->SetInput(smoother->GetOutput());
  writeDense->SetFileName(denseFileName);
  writeDense->Update();

  std::cout << "Shifting by " << threshold << std::endl;
  SubtractType::Pointer sub = SubtractType::New();
  sub->SetConstant1(threshold);
  sub->SetInput2(smoother->GetOutput());
  sub->Update();

  std::cout << "Writing out psi to " << psiFileName << std::endl;
  WriterType::Pointer writerPsi = WriterType::New();
  writerPsi->SetInput(sub->GetOutput());
  writerPsi->SetFileName(psiFileName);
  writerPsi->Update();

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

  std::cout << "Extending" << std::endl;
  SweepingType::Pointer sweeping = SweepingType::New();
  sweeping->SetNarrowbandInput(init->GetInitializedNarrowbandOutput());
  sweeping->SetMaskInput(init->GetNarrowbandMaskOutput());
  sweeping->SetMaxIteration(2);
  sweeping->SetMaxError(1e-11);
  sweeping->SetOrder(1);
  sweeping->Update();

  std::cout << "Writing out high order signed distance transform " << outputFileName << std::endl;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(sweeping->GetOutput());
  writer->SetFileName(outputFileName);
  writer->Update();

  return EXIT_SUCCESS;
}
