#include <iostream>
#include <map>
#include <fstream>
#include <experimental/filesystem>

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


using OutcomeMap = std::map< std::string, std::string >;
void write_to_csv(std::string filename, OutcomeMap outcomes) {
  std::string delim = ",";

  /* Does not exist, write header */
  if (!std::experimental::filesystem::exists(filename)) {
      std::ofstream file;
      file.open(filename, std::ofstream::out);
      if (file.is_open())
      {
        for (OutcomeMap::iterator it = outcomes.begin(); it != outcomes.end(); ++it) {
          file << it->first << delim;
        }
        file << std::endl;

        file.close();
      }
      else std::cerr << "Unable to open file " << filename << " for writing." << std::endl;
  }

  /* Write contents */
  std::ofstream file;
  file.open(filename, std::ofstream::out | std::ofstream::app);
  if (file.is_open())
  {
    for (OutcomeMap::iterator it = outcomes.begin(); it != outcomes.end(); ++it) {
      file << it->second << delim;
    }
    file << std::endl;

    file.close();
  }
  else std::cerr << "Unable to open file " << filename << " for appending." << std::endl;
}

int main(int argc, char * argv[])
{
  if( argc < 2 )
  {
    std::cerr << "Usage: "<< std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFileName> <CSVFileName>";
    std::cerr << std::endl;
    std::cerr << "The image is cropped and padded by `PadSize` to reduce computation time.";
    std::cerr << std::endl << std::endl;
    return EXIT_FAILURE;
  }

  /* Read input Parameters */
  std::string inputFileName = argv[1];
  std::string csvFileName = argv[2];

  std::cout << "Read in the following parameters:" << std::endl;
  std::cout << "  Input Filename: " << inputFileName << std::endl;
  std::cout << "  CSV Fileanme:   " << csvFileName << std::endl;
  std::cout << std::endl;

  /* Read */
  std::cout << "Reading in " << inputFileName << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputFileName);
  reader->Update();

  /* Mesh */
  std::cout << "Meshing" << std::endl;
  MesherType::Pointer meshSource = MesherType::New();
  meshSource->SetInput(reader->GetOutput());
  meshSource->Update();

  /* Volume */
  std::cout << "Computing" << std::endl;
  TConvert::Pointer convert = TConvert::New();
  convert->SetInput(meshSource->GetOutput());
  convert->Update();

  MeshCalcType::Pointer meshCalc = MeshCalcType::New();
  meshCalc->SetSimplexMesh(convert->GetOutput());
  meshCalc->Compute();
  std::cout << " Area:   " << meshCalc->GetArea() << std::endl;
  std::cout << " Volume: " << meshCalc->GetVolume() << std::endl;

  /* Outcomes */
  std::cout << "Writing to CSV file " << csvFileName << std::endl;
  OutcomeMap outcomes;
  outcomes["Algorithm"] = "Morphometry";
  outcomes["Input Filename"] = inputFileName;
  outcomes["CSV Filename"] = csvFileName;
  outcomes["Volume"] = std::to_string(meshCalc->GetVolume());
  outcomes["Area"] = std::to_string(meshCalc->GetArea());
  write_to_csv(csvFileName, outcomes);

  return EXIT_SUCCESS;
}
