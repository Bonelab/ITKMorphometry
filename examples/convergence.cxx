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

EmbeddingImageType::Pointer create_sphere(EmbeddingImageType::SpacingType physicalSize, EmbeddingImageType::SpacingType spacing, double radius)
{
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
  itk::ImageRegionIterator<EmbeddingImageType> in(image, region);
  for (in.GoToBegin(); !in.IsAtEnd(); ++in)
  {
    EmbeddingImageType::PointType thisIndex;
    image->TransformIndexToPhysicalPoint(in.GetIndex(), thisIndex);
    in.Set(thisIndex.EuclideanDistanceTo(center) - radius);
  }

  return image;
}

struct norm {
  double l1;
  double linf;
  double h;
};

norm measure_norm(EmbeddingImageType::Pointer const truth, EmbeddingImageType::Pointer const predicted, double h)
{
  /* Setup iterators */
  itk::ImageRegionConstIterator<EmbeddingImageType> t(truth, truth->GetRequestedRegion());
  itk::ImageRegionConstIterator<EmbeddingImageType> p(predicted, predicted->GetRequestedRegion());

  /* Compute best shift */
  int i = 0;
  double shift = 0;
  for (t.GoToBegin(), p.GoToBegin(); !t.IsAtEnd() && !p.IsAtEnd(); ++t, ++p)
  {
    double truth = t.Get();
    double predicted = p.Get();
    shift += (truth - predicted);
    i++;
  }
  shift /= (double)i;

  i = 0;
  double l1 = 0.;
  double linf = 0.;
  for (t.GoToBegin(), p.GoToBegin(); !t.IsAtEnd() && !p.IsAtEnd(); ++t, ++p)
  {
    double truth = t.Get();
    double predicted = p.Get() + shift;

    l1 += std::abs(truth - predicted);
    linf = std::max(linf, std::abs(truth - predicted));
    i++;
  }

  norm n;
  n.l1 = l1 / (double)i;
  n.linf = linf;
  n.h = h;
  return n;
}

norm run_value(double h, double radius)
{
  WriterType::Pointer writer = WriterType::New();

  EmbeddingImageType::SpacingType physicalSize(100.);
  // physicalSize[2] = h;
  EmbeddingImageType::SpacingType spacing(h);

  EmbeddingImageType::Pointer dt = create_sphere(physicalSize, spacing, radius);

  writer->SetFileName("dt.nii");
  writer->SetInput(dt);
  writer->Update();

  BinarizeType::Pointer thresh = BinarizeType::New();
  thresh->SetInput(dt);
  thresh->SetInsideValue(1);
  thresh->SetOutsideValue(0);
  thresh->SetUpperThreshold(0);

  Writer2Type::Pointer w2 = Writer2Type::New();
  w2->SetFileName("seg.nii");
  w2->SetInput(thresh->GetOutput());
  w2->Update();

  InitializeType::Pointer init = InitializeType::New();
  init->SetInput(thresh->GetOutput());
  init->SetLabel(1);
  init->SetDither(false);
  // init->Update();

  writer->SetFileName("init.nii");
  writer->SetInput(init->GetOutput());
  writer->Update();

  // return measure_norm(dt, init->GetOutput(), h);

  SDTType::Pointer correct = SDTType::New();
  // correct->SetInput(dt);
  correct->SetInput(init->GetOutput());
  correct->SetMaxIteration(100);
  correct->SetNarrowbandSize(15.);
  correct->SetMaxError(0.);
  correct->Update();

  writer->SetFileName("correct.nii");
  writer->SetInput(correct->GetOutput());
  writer->Update();

  return measure_norm(dt, correct->GetOutput(), h);
}

int main(int argc, char * argv[])
{
  std::string filename = "temp.nii";
  std::vector<double> hs {4., 2., 1.};//{8., 4., 2., 1., 0.5};// {2., 1., 0.5, 0.25};//{8., 4., 2., 1.};//0.25};//{2.0, 1.0, 0.5, 0.25};
  for(double &h: hs)
  {
    norm n = run_value(h, 25.0);
    std::cout << "! " << n.l1 << " " << n.linf << " " << n.h << std::endl;
  }

  return EXIT_SUCCESS;
}


// int main(int argc, char * argv[])
// {
//   std::string filename = "temp.nii";

//   EmbeddingImageType::SpacingType physicalSize(100.);
//   EmbeddingImageType::SpacingType spacing(1.);
//   double radius = 10.;

//   std::cout << "Creating sphere" << std::endl;
//   EmbeddingImageType::Pointer dt = create_sphere(physicalSize, spacing, radius);

//   std::cout << "Writing to " << filename << std::endl;
//   WriterType::Pointer writer = WriterType::New();
//   writer->SetFileName(filename);
//   writer->SetInput(dt);
//   writer->Update();

//   return EXIT_SUCCESS;
// }
