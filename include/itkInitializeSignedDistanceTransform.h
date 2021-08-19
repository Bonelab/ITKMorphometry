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
#ifndef itkInitializeSignedDistanceTransform_h
#define itkInitializeSignedDistanceTransform_h

#include "itkImageToImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkSignedDitherImageFilter.h"

namespace itk
{

/** \class InitializeSignedDistanceTransform
 *
 * \brief Initialize the signed distance transform to have no zeros
 *
 * \author Bryce Besler
 * \ingroup Morphometry
 */
template <typename TInputImage, typename TOutputImage>
class InitializeSignedDistanceTransform : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(InitializeSignedDistanceTransform);

  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;
  static constexpr unsigned int Order = 3;

  /** Type defines */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using InputPixelType = typename InputImageType::PixelType;
  using OutputPixelType = typename OutputImageType::PixelType;
  using RealType = typename NumericTraits< TOutputImage >::RealType;

  /** Standard class typedefs. */
  using Self = InitializeSignedDistanceTransform<InputImageType, TOutputImage>;
  using Superclass = ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(InitializeSignedDistanceTransform, ImageToImageFilter);

  /** Standard New macro. */
  itkNewMacro(Self);

  /** MaxNumberOfLabels macros */
  itkSetMacro(Label, InputPixelType);
  itkGetConstMacro(Label, InputPixelType);

  /** Set if should divide gamma. */
  itkSetMacro(Dither, bool);
  itkGetConstReferenceMacro(Dither, bool);

  /** Filter types */
  using ThresholdFilterType = BinaryThresholdImageFilter< InputImageType, InputImageType >;
  using SubtractFilterType = SubtractImageFilter< InputImageType, InputImageType, InputImageType >;
  using DistanceMapFilterType = SignedMaurerDistanceMapImageFilter< InputImageType, TOutputImage >;
  using SubtractSDTFilerType = SubtractImageFilter< TOutputImage, TOutputImage, TOutputImage >;
  using MultiplyFilterType = MultiplyImageFilter< TOutputImage, TOutputImage, TOutputImage >;
  using DitherFilterType = SignedDitherImageFilter< TOutputImage >;

protected:
  InitializeSignedDistanceTransform();
  ~InitializeSignedDistanceTransform() override = default;

  void GenerateData () override;

private:
  /* Member variables */
  bool m_Dither;
  InputPixelType m_Label;
  typename ThresholdFilterType::Pointer m_ThresholdFilter;
  typename SubtractFilterType::Pointer m_SubtractFilter;
  typename DistanceMapFilterType::Pointer m_ForegroundSDTFilter;
  typename DistanceMapFilterType::Pointer m_BackgroundSDTFilter;
  typename SubtractSDTFilerType::Pointer m_SubtractSDTFilter;
  typename MultiplyFilterType::Pointer m_MultiplyFilter;
  typename DitherFilterType::Pointer m_DitherFilter;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkInitializeSignedDistanceTransform.hxx"
#endif

#endif // itkInitializeSignedDistanceTransform_h
