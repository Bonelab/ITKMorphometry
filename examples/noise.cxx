#include <iostream>
#include <vector>
#include <string>

#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "itkInitializeSignedDistanceTransform.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkHighOrderSignedDistanceTransformCorrectionImageFilter.h"
#include "itkHighOrderSignedFastSweeping.h"

constexpr unsigned int ImageDimension = 3;
using EmbeddingPixelType = double;
using EmbeddingImageType = itk::Image<EmbeddingPixelType, ImageDimension>;
using BinaryPixelType = unsigned char;
using BinaryImageType = itk::Image<BinaryPixelType, ImageDimension>;

using BinarizeType = itk::BinaryThresholdImageFilter< EmbeddingImageType, BinaryImageType >;
using InitializeType = itk::InitializeSignedDistanceTransform< BinaryImageType, EmbeddingImageType >;
// using SDTType = itk::HighOrderSignedDistanceTransformCorrectionImageFilter< EmbeddingImageType >;
using SDTType = itk::HighOrderSignedFastSweeping< EmbeddingImageType >;

using WriterType = itk::ImageFileWriter<EmbeddingImageType>;
using Writer2Type = itk::ImageFileWriter<BinaryImageType>;

int main(int argc, char * argv[])
{
  double h = 0.1;
  double radius = 2.5;
  EmbeddingImageType::SpacingType physicalSize(10.);
  // physicalSize[2] = h;
  EmbeddingImageType::SpacingType spacing(h);
  int q = 3;

  /* Create region */
  EmbeddingImageType::IndexType start;
  EmbeddingImageType::SizeType size;
  EmbeddingImageType::PointType origin;
  EmbeddingImageType::PointType center;

  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    start[i] = 0;
    size[i] = int(physicalSize[i] / spacing[i]);
    origin[i] = -physicalSize[i] / 2.0;
    center[i] = -spacing[2];
  }

  EmbeddingImageType::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);

  /* Create image */
  typename EmbeddingImageType::Pointer image = EmbeddingImageType::New();
  image->SetRegions(region);
  image->SetOrigin(origin);
  image->SetSpacing(spacing);
  image->Allocate();

  /* Fill */
  double max = 0.;
  itk::ImageRegionIterator<EmbeddingImageType> in(image, region);
  for (in.GoToBegin(); !in.IsAtEnd(); ++in)
  {
    EmbeddingImageType::PointType thisIndex;
    image->TransformIndexToPhysicalPoint(in.GetIndex(), thisIndex);
    // std::cout << thisIndex << std::endl;
    // EmbeddingPixelType arg = 2*M_PI/h * (std::pow(thisIndex[0], 2) + std::pow(thisIndex[1], 2) + std::pow(thisIndex[2], 2));
    EmbeddingPixelType arg = 2*M_PI/(10.*h) * (thisIndex[0] + thisIndex[1] + thisIndex[2]) + h/2.;
    EmbeddingPixelType value = thisIndex.EuclideanDistanceTo(center) - radius + std::pow(h, q) * std::sin(arg);
    // std::cout << arg << " " << std::pow(h, q) * std::sin(arg) << " " << value << " " << thisIndex.EuclideanDistanceTo(center) - radius << std::endl;
    in.Set(value);
    // in.Set(std::pow(h, q) * std::sin(arg));
    max = std::max(max, std::abs(std::pow(h, q) * std::sin(arg)));
  }
  std::cout << "Max: " << max << std::endl;

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(image);
  writer->SetFileName("noise_q3.nii");
  writer->Update();

  return EXIT_SUCCESS;
}
