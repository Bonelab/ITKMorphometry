#include <iostream>
#include <vector>

// #include <random>

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

// distance = -1. * std::max(-1.*s1, -1.*s2);

      // /* Construct circle at intersection */
      // ImageType::PointType::VectorType z = center2.GetVectorFromOrigin() - center1.GetVectorFromOrigin();
      // double d = z.Normalize();
      // double r = (-d + radius - radius) * (-d - radius + radius) * (-d + radius + radius) * (d + radius + radius);
      // r = std::sqrt(r) / (2.*d);

      // /* Distance to circle */
      // ImageType::PointType::VectorType c = 0.5 * (center1.GetVectorFromOrigin() + center2.GetVectorFromOrigin());
      // ImageType::PointType::VectorType x = thisIndex.GetVectorFromOrigin() - c;
      // ImageType::PointType::VectorType y = x - (x * z)*z;
      // y.Normalize();
      // ImageType::PointType::VectorType q = c + r*y;

      // // double a = r - q.GetNorm();
      // // double b = (x - y).GetNorm();
      // distance = -1. * (q - x).GetNorm();

      // /* Construct circle at intersection */
      // ImageType::PointType::VectorType z = center2.GetVectorFromOrigin() - center1.GetVectorFromOrigin();
      // double d = z.Normalize();
      // double r = std::sqrt(radius*radius - (d/2.)*(d/2.));
      // // double r = (-d + radius - radius) * (-d - radius + radius) * (-d + radius + radius) * (d + radius + radius);
      // // double r = d*d - radius*radius + radius*radius;
      // // r = 4.*d*d*radius*radius - r*r;
      // // r = std::sqrt(r) / (2.*d);

      // /* Distance to circle */
      // ImageType::PointType::VectorType c = 0.5 * (center1.GetVectorFromOrigin() + center2.GetVectorFromOrigin());
      // ImageType::PointType::VectorType x = thisIndex.GetVectorFromOrigin();
      // ImageType::PointType::VectorType t = x - c;
      // ImageType::PointType::VectorType y = t - (t * z)*z;

      // double a = r - y.GetNorm();
      // double b = (x - y).GetNorm();
      // distance = -1. * std::sqrt(a*a + b*b);



      // /* Construct circle at intersection */
      // ImageType::PointType::VectorType n = center2.GetVectorFromOrigin() - center1.GetVectorFromOrigin();
      // double d = n.Normalize();
      // double r = std::sqrt(radius*radius - (d/2.)*(d/2.))/2.;
      // ImageType::PointType::VectorType c = 0.5 * (center1.GetVectorFromOrigin() + center2.GetVectorFromOrigin());
      // ImageType::PointType::VectorType p = thisIndex.GetVectorFromOrigin() - c;

      // double a = r - CrossProduct(n,p).GetNorm();
      // double b = p*n;
      // distance = -1. * std::sqrt(a*a + b*b);

  return image;
}

struct norm {
  double l1;
  double linf;
};

norm measure_norm(ImageType::Pointer truth, ImageType::Pointer predicted, BinaryImageType::Pointer mask, double h)
{
  /* Setup iterators */
  itk::ImageRegionConstIterator<ImageType> t(truth, truth->GetRequestedRegion());
  itk::ImageRegionConstIterator<ImageType> p(predicted, predicted->GetRequestedRegion());
  itk::ImageRegionConstIterator<BinaryImageType> m(mask, mask->GetRequestedRegion());

  // std::default_random_engine generator;
  // std::uniform_real_distribution<double> distribution(-1.0,1.0);

  /* Perform */
  unsigned int i = 0;
  double l1 = 0.;
  double linf = 0.;
  for (
    t.GoToBegin(), p.GoToBegin(), m.GoToBegin();
    !t.IsAtEnd() && !p.IsAtEnd() && !m.IsAtEnd();
    ++t, ++p, ++m
  )
  {
    if (m.Get() == 0) {
      continue;
    }

    double truth = t.Get();
    double predicted = p.Get();
    // double predicted = truth + h*h * distribution(generator);

    l1 += std::abs(truth - predicted);
    linf = std::max(linf, std::abs(truth - predicted));
    i++;
  }

  norm n;
  n.l1 = l1 / (double)i;
  n.linf = linf;
  return n;
}

norm compute_error(ImageType::Pointer dt, double rho1, double rho2, double epsilon, double h) {
  char filename4[100];
  sprintf(filename4, "narrow_band_dt_%f.nii", h);
  WriterType::Pointer writer4 = WriterType::New();
  writer4->SetFileName(filename4);
  writer4->SetInput(dt);
  writer4->Update();

  // Density image
  CombineFilter::Pointer combine = CombineFilter::New();
  combine->SetInput(dt);
  combine->SetRho1(rho1);
  combine->SetRho2(rho2);
  combine->SetEpsilon(epsilon);

  // Shift
  double threshold = (rho1 + rho2) / 2.;
  SubtractType::Pointer sub = SubtractType::New();
  sub->SetConstant1(threshold);
  sub->SetInput2(combine->GetOutput());

  char filename3[100];
  sprintf(filename3, "narrow_band_density_%f.nii", h);
  WriterType::Pointer writer3 = WriterType::New();
  writer3->SetFileName(filename3);
  writer3->SetInput(sub->GetOutput());
  writer3->Update();

  // NB select
  NBType::Pointer nbSelect = NBType::New();
  nbSelect->SetInput(sub->GetOutput());
  // nbSelect->SetInput(dt);
  nbSelect->SetRadius(3);
  nbSelect->Update();

  // Init
  InitType::Pointer init = InitType::New();
  init->SetGrayInput(sub->GetOutput());
  // init->SetGrayInput(dt);
  init->SetNarrowbandInput(nbSelect->GetOutput());
  init->SetMaxIterations(100);
  init->SetConvergenceOrder(3);
  init->Update();

  char filename[100];
  sprintf(filename, "narrow_band_solved_%f.nii", h);
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(init->GetInitializedNarrowbandOutput());
  writer->Update();

  sprintf(filename, "narrow_band_mask_%f.nii", h);
  BWriterType::Pointer writer2 = BWriterType::New();
  writer2->SetFileName(filename);
  writer2->SetInput(init->GetNarrowbandMaskOutput());
  writer2->Update();

  return measure_norm(dt, init->GetInitializedNarrowbandOutput(), init->GetNarrowbandMaskOutput(), h);
}

int main() {
  double r2 = 5.;
  double c1 = 1.75;
  double c2 = -1.75;
  double epsilon = 2.;
  double rho1 = 0.;
  double rho2 = 100.;
  double physicalSize = 20.;
  double radius = 5.;
  double a = 2.;
  double c = 4.;
  // std::vector<float> hVector{0.5};
  std::vector<float> hVector{0.5, 0.25, 0.125, 0.0625};
  // std::vector<float> hVector{0.2, 0.1, 0.05, 0.025};

  // std::cout << "Sphere" << std::endl;
  // for(auto &h: hVector) {
  //   // Generate dt
  //   ImageType::SpacingType physicalSizeArray(physicalSize);
  //   ImageType::SpacingType spacing(h);
  //   ImageType::Pointer dt = create_sphere(physicalSizeArray, spacing, radius);

  //   // Measure norm
  //   norm n = compute_error(dt, rho1, rho2, epsilon, h);

  //   // char filename[100];
  //   // sprintf(filename, "narrow_band_%f.nii", h);
  //   // WriterType::Pointer writer = WriterType::New();
  //   // writer->SetFileName(filename);
  //   // writer->SetInput(dt);
  //   // writer->Update();

  //   std::cout << h << " " << n.l1 << " " << n.linf << std::endl << std::flush;
  // }

  // std::cout << "Torus" << std::endl;
  // for(auto &h: hVector) {
  //   // Generate dt
  //   ImageType::SpacingType physicalSizeArray(physicalSize);
  //   ImageType::SpacingType spacing(h);
  //   ImageType::Pointer dt = create_torus(physicalSizeArray, spacing, a, c);

  //   // Measure norm
  //   norm n = compute_error(dt, rho1, rho2, epsilon, h);

  //   // char filename[100];
  //   // sprintf(filename, "narrow_band_torus_%f.nii", h);
  //   // WriterType::Pointer writer = WriterType::New();
  //   // writer->SetFileName(filename);
  //   // writer->SetInput(dt);
  //   // writer->Update();

  //   std::cout << h << " " << n.l1 << " " << n.linf << std::endl << std::flush;
  // }


  // std::cout << "Spindle Torus" << std::endl;
  // for(auto &h: hVector) {
  //   // Generate dt
  //   ImageType::SpacingType physicalSizeArray(physicalSize);
  //   ImageType::SpacingType spacing(h);
  //   ImageType::Pointer dt = create_torus(physicalSizeArray, spacing, 4., 3.5);

  //   // Measure norm
  //   norm n = compute_error(dt, rho1, rho2, epsilon, h);

  //   // char filename[100];
  //   // sprintf(filename, "narrow_band_torus_%f.nii", h);
  //   // WriterType::Pointer writer = WriterType::New();
  //   // writer->SetFileName(filename);
  //   // writer->SetInput(dt);
  //   // writer->Update();

  //   std::cout << h << " " << n.l1 << " " << n.linf << std::endl << std::flush;
  // }

  std::cout << "Double sphere" << std::endl;
  for(auto &h: hVector) {
    // Generate dt
    ImageType::SpacingType physicalSizeArray(physicalSize);
    ImageType::SpacingType spacing(h);
    ImageType::Pointer dt = create_double_sphere(physicalSizeArray, spacing, c1, c2, r2);

    // Measure norm
    norm n = compute_error(dt, rho1, rho2, epsilon, h);

    // char filename[100];
    // sprintf(filename, "narrow_band_ds_%f.nii", h);
    // WriterType::Pointer writer = WriterType::New();
    // writer->SetFileName(filename);
    // writer->SetInput(dt);
    // writer->Update();

    std::cout << h << " " << n.l1 << " " << n.linf << std::endl << std::flush;
  }

  return EXIT_SUCCESS;
}
