#include <iostream>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCleanBinarySurfaceImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkSimplexMeshVolumeCalculator.h"
#include "itkTriangleMeshToSimplexMeshFilter.h"

/* Setup Types */
constexpr unsigned int ImageDimension = 3;
using InputPixelType = unsigned int;
using InputImageType = itk::Image<InputPixelType, ImageDimension>;
using OutputPixelType = InputPixelType;
using OutputImageType = itk::Image<OutputPixelType, ImageDimension>;

using ReaderType = itk::ImageFileReader< InputImageType >;
using WriterType = itk::ImageFileWriter< OutputImageType >;
using CleanBinarySurfaceImageFilterType = itk::CleanBinarySurfaceImageFilter< OutputImageType >;
using BinaryThresholdImageFilterType = itk::BinaryThresholdImageFilter< InputImageType, OutputImageType >;
using MeshType = itk::Mesh<float, 3>;
using MesherType = itk::BinaryMask3DMeshSource< OutputImageType, MeshType >;
using TSimplex = itk::SimplexMesh<float, 3>;
using TConvert = itk::TriangleMeshToSimplexMeshFilter<MeshType, TSimplex>;
using MeshCalcType = itk::SimplexMeshVolumeCalculator< TSimplex >;

int main(int argc, char * argv[])
{
  if( argc < 4 )
  {
    std::cerr << "Usage: "<< std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFileName> <OutputMeasure> <Label> <PadSize>";
    std::cerr << std::endl;
    std::cerr << "The image is cropped and padded by `PadSize` to reduce computation time.";
    std::cerr << std::endl << std::endl;
    return EXIT_FAILURE;
  }

  /* Read input Parameters */
  std::string inputFileName = argv[1];
  std::string outputFileName = argv[2];
  int label = std::stoi(argv[3]);
  int padSize = std::stoi(argv[4]);

  std::cout << "Read in the following parameters:" << std::endl;
  std::cout << "  Input Filename:  " << inputFileName << std::endl;
  std::cout << "  Output Fileanme: " << outputFileName << std::endl;
  std::cout << "  Label:           " << label << std::endl;
  std::cout << "  Pad size:        " << padSize << std::endl;
  std::cout << std::endl;

  /* Read */
  std::cout << "Reading in " << inputFileName << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputFileName);
  reader->Update();

  /* Threshold */
  std::cout << "Thresholding at " << label << std::endl;
  BinaryThresholdImageFilterType::Pointer thresh = BinaryThresholdImageFilterType::New();
  thresh->SetInput(reader->GetOutput());
  thresh->SetLowerThreshold(label);
  thresh->SetUpperThreshold(label);
  thresh->SetInsideValue(1);
  thresh->SetOutsideValue(0);
  thresh->Update();

  /* Mesh */
  std::cout << "Meshing" << std::endl;
  MesherType::Pointer meshSource = MesherType::New();
  meshSource->SetInput(thresh->GetOutput());
  meshSource->Update();

  /* Volume */
  std::cout << "Computing" << std::endl;
  TConvert::Pointer convert = TConvert::New();
  convert->SetInput(meshSource->GetOutput());
  convert->Update();

  MeshCalcType::Pointer meshCalc = MeshCalcType::New();
  meshCalc->SetSimplexMesh(convert->GetOutput());
  meshCalc->Compute();
  std::cout << " Area:" << meshCalc->GetArea() << std::endl;
  std::cout << " Volume:" << meshCalc->GetVolume() << std::endl;

  return EXIT_SUCCESS;
}
