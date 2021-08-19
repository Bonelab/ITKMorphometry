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
#ifndef itkSmoothBinarySurfaceImageFilter_h
#define itkSmoothBinarySurfaceImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{

/** \class SmoothBinarySurfaceImageFilter
 *
 * \brief High order signed fast sweeping method.
 *
 * \author Bryce Besler
 * \ingroup Morphometry
 */
template <typename TInputImage>
class SmoothBinarySurfaceImageFilter : public ImageToImageFilter<TInputImage, TInputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(SmoothBinarySurfaceImageFilter);

  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int Order = 3;

  /** Type defines */
  using InputImageType = TInputImage;
  using InputPixelType = typename InputImageType::PixelType;
  using RealType = typename NumericTraits< InputPixelType >::RealType;
  using OffsetType = typename InputImageType::OffsetType;
  using IndexType = typename TInputImage::IndexType;
  using StencilType = itk::Vector< OffsetType, 2*Order >;
  using StencilsType = itk::Vector< StencilType, InputImageDimension >;
  using ImageRegionType = typename TInputImage::RegionType;

  /** Standard class typedefs. */
  using Self = SmoothBinarySurfaceImageFilter<InputImageType>;
  using Superclass = ImageToImageFilter<InputImageType, TInputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(SmoothBinarySurfaceImageFilter, ImageToImageFilter);

  /** Standard New macro. */
  itkNewMacro(Self);

protected:
  SmoothBinarySurfaceImageFilter();
  ~SmoothBinarySurfaceImageFilter() override = default;

  // Override since the filter needs all the data for the algorithm.
  void GenerateInputRequestedRegion() override;

  // Override since the filter produces the entire dataset.
  void EnlargeOutputRequestedRegion(DataObject * output) override;

  RealType SolveQuadratic(IndexType &index);
  StencilsType GenerateStencils();

  // void BeforeThreadedGenerateData() override;
  void AfterThreadedGenerateData() override;
  void DynamicThreadedGenerateData (const ImageRegionType &outputRegionForThread) override;

private:

};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSmoothBinarySurfaceImageFilter.hxx"
#endif

#endif // itkSmoothBinarySurfaceImageFilter_h
