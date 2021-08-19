#include <iostream>
#include <vector>

#include "itkImageFileWriter.h"
#include "itkInitializeSignedDistanceTransform.h"
#include "itkCombineTwoPhaseImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkSelectNarrowbandImageFilter.h"
#include "itkInitializeTwoPhaseNarrowBandImageFilter.h"

constexpr unsigned int ImageDimension = 3;
using PixelType = double;
using ImageType = itk::Image<PixelType, ImageDimension>;
using BinaryPixelType = unsigned int;
using BinaryImageType = itk::Image<BinaryPixelType, ImageDimension>;

using BWriterType = itk::ImageFileWriter<BinaryImageType>;
using WriterType = itk::ImageFileWriter<ImageType>;
using CombineFilter = itk::CombineTwoPhaseImageFilter<ImageType>;
using SubtractType = itk::SubtractImageFilter< ImageType, ImageType, ImageType >;
using NBType = itk::SelectNarrowbandImageFilter<ImageType, BinaryImageType>;
using InitType = itk::InitializeTwoPhaseNarrowBandImageFilter<ImageType, BinaryImageType>;

ImageType::Pointer create_sphere(ImageType::SpacingType physicalSize, ImageType::SpacingType spacing, double radius)
{
  /* Create region */
  ImageType::IndexType start;
  ImageType::SizeType size;
  ImageType::PointType origin;
  ImageType::PointType center;

  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    start[i] = 0;
    size[i] = int(physicalSize[i] / spacing[i]);
    origin[i] = -physicalSize[i] / 2.0;
    center[i] = -spacing[2];
  }

  ImageType::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);

  /* Create image */
  typename ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->SetOrigin(origin);
  image->SetSpacing(spacing);
  image->Allocate();

  /* Fill */
  itk::ImageRegionIterator<ImageType> in(image, region);
  for (in.GoToBegin(); !in.IsAtEnd(); ++in)
  {
    ImageType::PointType thisIndex;
    image->TransformIndexToPhysicalPoint(in.GetIndex(), thisIndex);
    in.Set(thisIndex.EuclideanDistanceTo(center) - radius);
  }

  return image;
}

ImageType::Pointer create_torus(ImageType::SpacingType physicalSize, ImageType::SpacingType spacing, double a, double c)
{
  /* Create region */
  ImageType::IndexType start;
  ImageType::SizeType size;
  ImageType::PointType origin;
  ImageType::PointType center;

  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    start[i] = 0;
    size[i] = int(physicalSize[i] / spacing[i]);
    origin[i] = -physicalSize[i] / 2.0;
    center[i] = -spacing[2];
  }

  ImageType::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);

  /* Create image */
  typename ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->SetOrigin(origin);
  image->SetSpacing(spacing);
  image->Allocate();

  /* Fill */
  itk::ImageRegionIterator<ImageType> in(image, region);
  for (in.GoToBegin(); !in.IsAtEnd(); ++in)
  {
    ImageType::PointType thisIndex;
    image->TransformIndexToPhysicalPoint(in.GetIndex(), thisIndex);
    ImageType::PointType d = thisIndex - center;
    
    double q = std::sqrt(d[0]*d[0] + d[1]*d[1]);
    q = q - c;
    q = q*q + d[2]*d[2];
    in.Set(std::sqrt(q) - a);
  }

  return image;
}

ImageType::Pointer create_double_sphere(ImageType::SpacingType physicalSize, ImageType::SpacingType spacing, double c1, double c2, double radius)
{
  /* Create region */
  ImageType::IndexType start;
  ImageType::SizeType size;
  ImageType::PointType origin;
  ImageType::PointType center1;
  ImageType::PointType center2;

  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    start[i] = 0;
    size[i] = int(physicalSize[i] / spacing[i]);
    origin[i] = -physicalSize[i] / 2.0;
    center1[i] = c1;
    center2[i] = c2;
  }

  ImageType::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);

  /* Create image */
  typename ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->SetOrigin(origin);
  image->SetSpacing(spacing);
  image->Allocate();

  /* Fill */
  itk::ImageRegionIterator<ImageType> in(image, region);
  for (in.GoToBegin(); !in.IsAtEnd(); ++in)
  {
    ImageType::PointType thisIndex;
    image->TransformIndexToPhysicalPoint(in.GetIndex(), thisIndex);

    /* Sample spheres*/
    double d1 = thisIndex.EuclideanDistanceTo(center1);
    double d2 = thisIndex.EuclideanDistanceTo(center2);
    double outside1 = d1 > radius;
    double outside2 = d2 > radius;
    double s1 = thisIndex.EuclideanDistanceTo(center1) - radius;
    double s2 = thisIndex.EuclideanDistanceTo(center2) - radius;
  
    /* Intersecting sphere */
    double intersectionDistance = 0.;
    ImageType::PointType::VectorType n = center2.GetVectorFromOrigin() - center1.GetVectorFromOrigin();
    double d = n.Normalize();
    double r = std::sqrt(radius*radius - (d/2.)*(d/2.));

    ImageType::PointType::VectorType c = 0.5 * (center1.GetVectorFromOrigin() + center2.GetVectorFromOrigin());
    ImageType::PointType::VectorType x = thisIndex.GetVectorFromOrigin();
    ImageType::PointType::VectorType t = x - c;
    ImageType::PointType::VectorType y = t - (t * n)*n;

    double a = r - y.GetNorm();
    double b = t*n;
    intersectionDistance = -1. * std::sqrt(a*a + b*b);

    /* Test if inside cone */
    double tt = d/2. - std::abs(t*n);
    double hh = y.GetNorm();
    double h = tt * r / (d/2.);
    bool in_center = true;
    if (hh > h) {
      in_center = false;
    }

    /* Compute distance */
    double distance = 0.;
    if (outside1 && outside2) {
      distance = std::min(s1, s2);
    } else if (outside1 && !outside2) {
      if (in_center) {
        distance = intersectionDistance;
      } else {
        distance = s2;
      }
    } else if (!outside1 && outside2) {
      if (in_center) {
        distance = intersectionDistance;
      } else {
        distance = s1;
      }
    } else {
      distance = intersectionDistance;
    }
    in.Set(distance);
  }

  return image;
}

int main() {
  double r2 = 5.;
  double c1 = 1.75;
  double c2 = -1.75;
  double epsilon = 2.;
  double rho1 = 0.;
  double rho2 = 100.;
  double physicalSize = 30.;
  double radius = 5.;
  double a = 2.;
  double c = 4.;

  double min_ratio = 0.5;
  double max_ratio = 5.0;
  int number = 100;
  double slope = (max_ratio - min_ratio) / (number - 1);
  double inter = min_ratio;
  char filename[100];

  std::cout << "Sphere" << std::endl;
  for (unsigned int u = 0; u < number; u++) {
    double this_ratio = slope*u + inter;
    double h = radius / this_ratio;

    // Generate dt
    ImageType::SpacingType physicalSizeArray(physicalSize);
    ImageType::SpacingType spacing(h);
    ImageType::Pointer dt = create_sphere(physicalSizeArray, spacing, radius);

    sprintf(filename, "./rho_dir/dt_sphere_%f.nii", this_ratio);
    WriterType::Pointer writer2 = WriterType::New();
    writer2->SetFileName(filename);
    writer2->SetInput(dt);
    writer2->Update();

    // Density image
    CombineFilter::Pointer combine = CombineFilter::New();
    combine->SetInput(dt);
    combine->SetRho1(rho1);
    combine->SetRho2(rho2);
    combine->SetEpsilon(epsilon);

    sprintf(filename, "./rho_dir/density_sphere_%f.nii", this_ratio);
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(filename);
    writer->SetInput(combine->GetOutput());
    writer->Update();
  }

  std::cout << "Torus" << std::endl;
  for (unsigned int u = 0; u < number; u++) {
    double this_ratio = slope*u + inter;
    double h = a / this_ratio;

    // Generate dt
    ImageType::SpacingType physicalSizeArray(physicalSize);
    ImageType::SpacingType spacing(h);
    ImageType::Pointer dt = create_torus(physicalSizeArray, spacing, a, c);

    sprintf(filename, "./rho_dir/dt_torus_%f.nii", this_ratio);
    WriterType::Pointer writer2 = WriterType::New();
    writer2->SetFileName(filename);
    writer2->SetInput(dt);
    writer2->Update();

    // Density image
    CombineFilter::Pointer combine = CombineFilter::New();
    combine->SetInput(dt);
    combine->SetRho1(rho1);
    combine->SetRho2(rho2);
    combine->SetEpsilon(epsilon);

    char filename[100];
    sprintf(filename, "./rho_dir/density_torus_%f.nii", this_ratio);
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(filename);
    writer->SetInput(combine->GetOutput());
    writer->Update();
  }

  std::cout << "Double sphere" << std::endl;
  for (unsigned int u = 0; u < number; u++) {
    double this_ratio = slope*u + inter;
    double h = r2 / this_ratio;

    // Generate dt
    ImageType::SpacingType physicalSizeArray(physicalSize);
    ImageType::SpacingType spacing(h);
    ImageType::Pointer dt = create_double_sphere(physicalSizeArray, spacing, c1, c2, r2);

    sprintf(filename, "./rho_dir/dt_double_sphere_%f.nii", this_ratio);
    WriterType::Pointer writer2 = WriterType::New();
    writer2->SetFileName(filename);
    writer2->SetInput(dt);
    writer2->Update();

    // Density image
    CombineFilter::Pointer combine = CombineFilter::New();
    combine->SetInput(dt);
    combine->SetRho1(rho1);
    combine->SetRho2(rho2);
    combine->SetEpsilon(epsilon);

    char filename[100];
    sprintf(filename, "./rho_dir/density_double_sphere_%f.nii", this_ratio);
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(filename);
    writer->SetInput(combine->GetOutput());
    writer->Update();
  }

  return EXIT_SUCCESS;
}
