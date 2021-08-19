/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkInitializeTwoPhaseNarrowBandImageFilter_h
#define itkInitializeTwoPhaseNarrowBandImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkSelectNarrowbandImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"

namespace itk
{

/** \class InitializeTwoPhaseNarrowBandImageFilter
 *
 * \brief Compute a signed dither of an image.
 *
 * The filter guaranttes the sign of the embedding does not change
 * when dithering.
 * 
 *  Epsilon           Computed convergence criterion
 *  ConvergenceOrder  What power of spacing is used for convergence
 *  SplineOrder       BSpline order (0 - 5)
 *  MaxIterations     Maximum number of iterations
 *
 * Epsilon = min(spacing) ^ ConvergenceOrder
 * 
 * \author Bryce Besler
 * \ingroup Morphometry
 */
template <typename TGrayImage, typename TNarrowbandImage>
class InitializeTwoPhaseNarrowBandImageFilter : public ImageToImageFilter<TGrayImage, TGrayImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(InitializeTwoPhaseNarrowBandImageFilter);

  /** ImageDimension constants */
  itkStaticConstMacro(GrayImageDimension, SizeValueType, TGrayImage::ImageDimension);
  itkStaticConstMacro(NarrowbandImageDimension, SizeValueType, TNarrowbandImage::ImageDimension);

  /** Grey type defines */
  using GrayPixelType = typename TGrayImage::PixelType;
  using RegionType = typename TGrayImage::RegionType;
  using RealType = typename NumericTraits< GrayPixelType >::RealType;

  /** Narrowband type defines */
  using NarrowbandPixelType = typename TNarrowbandImage::PixelType;

  /** Standard class typedefs. */
  using Self = InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>;
  using Superclass = ImageToImageFilter<TGrayImage, TGrayImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(InitializeTwoPhaseNarrowBandImageFilter, ImageToImageFilter);

  /** Standard New macro. */
  itkNewMacro(Self);

  /** Image setters */
  void SetGrayInput(const TGrayImage * image1);
  const TGrayImage * GetGrayInput();
  void SetNarrowbandInput(const TNarrowbandImage * image1);
  const TNarrowbandImage * GetNarrowbandInput();

  /** Allow multiple outputs */
  using DataObjectPointer = typename Superclass::DataObjectPointer;
  using DataObjectPointerArraySizeType = ProcessObject::DataObjectPointerArraySizeType;
  // using Superclass::MakeOutput;
  DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx) override;

  TGrayImage * GetInitializedNarrowbandOutput();
  TNarrowbandImage * GetNarrowbandMaskOutput();

  /** Interpolator */
  using InterpolateType = BSplineInterpolateImageFunction< TGrayImage, RealType, RealType >;
  using PointType = typename InterpolateType::PointType;
  using VectorType = Vector< RealType, GrayImageDimension>;

  /** Epsilon macros */
  itkGetConstMacro(Epsilon, RealType);

  /** SplineOrder macros */
  itkSetMacro(SplineOrder, SizeValueType);
  itkGetConstMacro(SplineOrder, SizeValueType);

  /** ConvergenceOrder macros */
  itkSetMacro(ConvergenceOrder, SizeValueType);
  itkGetConstMacro(ConvergenceOrder, SizeValueType);

  /** MaxIterations macros */
  itkSetMacro(MaxIterations, SizeValueType);
  itkGetConstMacro(MaxIterations, SizeValueType);

#ifdef ITK_USE_CONCEPT_CHECKING
  // Begin concept checking
  itkConceptMacro(SameDimensionCheck,
                  (Concept::SameDimension<itkGetStaticConstMacro(GrayImageDimension),
                                          itkGetStaticConstMacro(NarrowbandImageDimension)>));
  // End concept checking
#endif

protected:
  InitializeTwoPhaseNarrowBandImageFilter();
  ~InitializeTwoPhaseNarrowBandImageFilter() override = default;

  // void PrintSelf(std::ostream & os, Indent indent) const override;

  void BeforeThreadedGenerateData() override;
  void DynamicThreadedGenerateData(const RegionType & region) override;
  void AfterThreadedGenerateData() override;

  void FindClosestPoint(PointType &inputPoint, PointType &closestPoint, SizeValueType& iterations);
  void FindClosestPointFirstOrder(PointType &point, SizeValueType& iterations);
  void FindClosestPointPerpendicular(PointType &point, SizeValueType& iterations);
  VectorType VelocityUpdate(const PointType &point, GrayPixelType &value, GrayPixelType &norm);

private:
  /* Member variables */
  RealType m_Epsilon;
  SizeValueType m_SplineOrder;
  SizeValueType m_ConvergenceOrder;
  SizeValueType m_MaxIterations;
  typename InterpolateType::Pointer m_Interpolator;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkInitializeTwoPhaseNarrowBandImageFilter.hxx"
#endif

#endif // itkInitializeTwoPhaseNarrowBandImageFilter_h
