#include <iostream>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCombineTwoPhaseImageFilter.h"

/* Setup Types */
constexpr unsigned int ImageDimension = 3;
using PixelType = double;
using ImageType = itk::Image<PixelType, ImageDimension>;

using ReaderType = itk::ImageFileReader< ImageType >;
using WriterType = itk::ImageFileWriter< ImageType >;
using CombineType = itk::CombineTwoPhaseImageFilter< ImageType >;

int main(int argc, char * argv[])
{
  if( argc < 7 )
  {
    std::cerr << "Usage: "<< std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFileName> <OutputFileName> <Epsilon> <Rho1> <Rho2> <Sin|Tanh>";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  /* Read input Parameters */
  std::string inputFileName = argv[1];
  std::string outputFileName = argv[2];
  double epsilon = std::stof(argv[3]);
  double rho1 = std::stof(argv[4]);
  double rho2 = std::stof(argv[5]);
  std::string type = argv[6];

  std::cout << "Read in the following parameters:" << std::endl;
  std::cout << "  Input Filename:  " << inputFileName << std::endl;
  std::cout << "  Output Fileanme: " << outputFileName << std::endl;
  std::cout << "  Epsilon:         " << epsilon << std::endl;
  std::cout << "  Rho1:            " << rho1 << std::endl;
  std::cout << "  Rho2:            " << rho2 << std::endl;
  std::cout << "  Type:            " << type << std::endl;
  std::cout << std::endl;

  /* Read */
  std::cout << "Reading in " << inputFileName << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputFileName);
  reader->Update();

  /* Combining */
  std::cout << "Combining" << std::endl;
  CombineType::Pointer combiner = CombineType::New();
  combiner->SetInput(reader->GetOutput());
  combiner->SetRho1(rho1);
  combiner->SetRho2(rho2);
  combiner->SetEpsilon(epsilon);
  if (type == "Sin") {
    combiner->SetApproximationToSin();
  } else if (type == "Tanh") {
    combiner->SetApproximationToTanh();
  } else {
    std::cerr << "Bad type: " << type << std::endl;
    return EXIT_FAILURE;
  }
  combiner->Update();

  /* Write */
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(combiner->GetOutput());
  writer->SetFileName(outputFileName);

  std::cout << "Writing results to " << outputFileName << std::endl;
  writer->Write();

  return EXIT_SUCCESS;
}
