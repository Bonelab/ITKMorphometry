#include <iostream>
#include <chrono>
#include <map>
#include <fstream>
#include <experimental/filesystem>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"

/* Setup Types */
constexpr unsigned int ImageDimension = 3;
using InputPixelType = unsigned int;
using InputImageType = itk::Image<InputPixelType, ImageDimension>;

using ReaderType = itk::ImageFileReader< InputImageType >;
using WriterType = itk::ImageFileWriter< InputImageType >;
using StructuringElementType = itk::BinaryBallStructuringElement< InputPixelType, ImageDimension >;
using ClosingType = itk::BinaryMorphologicalClosingImageFilter< InputImageType, InputImageType, StructuringElementType >;
using OpenType = itk::BinaryMorphologicalOpeningImageFilter< InputImageType, InputImageType, StructuringElementType >;

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
  if( argc < 4 )
  {
    std::cerr << "Usage: "<< std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFileName> <OutputFileName> <CSVFileName>";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  /* Read input Parameters */
  std::string inputFileName = argv[1];
  std::string outputFileName = argv[2];
  std::string csvFileName = argv[3];

  std::cout << "Read in the following parameters:" << std::endl;
  std::cout << "  Input Filename:  " << inputFileName << std::endl;
  std::cout << "  Output Filename: " << outputFileName << std::endl;
  std::cout << "  CSV Filename:    " << csvFileName << std::endl;
  std::cout << std::endl;

  /* Read */
  std::cout << "Reading in " << inputFileName << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputFileName);
  reader->Update();

  /* Clean */
  std::cout << "Morphological closing!" << std::endl;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(1);
  structuringElement.CreateStructuringElement();

  ClosingType::Pointer closer = ClosingType::New();
  closer->SetInput(reader->GetOutput());
  closer->SetForegroundValue(1);
  closer->SetKernel(structuringElement);

  OpenType::Pointer opener = OpenType::New();
  opener->SetInput(closer->GetOutput());
  opener->SetForegroundValue(1);
  opener->SetKernel(structuringElement);

  auto start = std::chrono::high_resolution_clock::now(); 
  opener->Update();
  auto stop = std::chrono::high_resolution_clock::now(); 
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "Exectuion time [us]: " << duration.count() << std::endl; 

  /* Write */
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(opener->GetOutput());
  writer->SetFileName(outputFileName);

  std::cout << "Writing results to " << outputFileName << std::endl;
  writer->Write();

  /* Outcomes */
  OutcomeMap outcomes;
  outcomes["Algorithm"] = "Open(Close(I))";
  outcomes["Input Filename"] = inputFileName;
  outcomes["Output Filename"] = outputFileName;
  outcomes["CSV Filename"] = csvFileName;
  outcomes["Execution time [us]"] = std::to_string(duration.count());
  std::cout << "Writing to CSV file " << csvFileName << std::endl;
  write_to_csv(csvFileName, outcomes);

  return EXIT_SUCCESS;
}
