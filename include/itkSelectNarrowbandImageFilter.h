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
#ifndef itkSelectNarrowbandImageFilter_h
#define itkSelectNarrowbandImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryCrossStructuringElement.h"

namespace itk
{

/** \class SelectNarrowbandImageFilter
 *
 * \brief Morphological dilation and erosion to select the narrowband voxels.
 *
 * \author Bryce Besler
 * \ingroup Morphometry
 */
template <typename TInputImage, typename TOutputImage>
class SelectNarrowbandImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(SelectNarrowbandImageFilter);

  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

  /** Type defines */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using InputPixelType = typename InputImageType::PixelType;
  using OutputPixelType = typename OutputImageType::PixelType;

  /** Standard class typedefs. */
  using Self = SelectNarrowbandImageFilter<InputImageType, TOutputImage>;
  using Superclass = ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(SelectNarrowbandImageFilter, ImageToImageFilter);

  /** Standard New macro. */
  itkNewMacro(Self);

  /** Radius macros */
  itkSetMacro(Radius, SizeValueType);
  itkGetConstMacro(Radius, SizeValueType);

  /** Filter types */
  using TKernel = BinaryCrossStructuringElement<OutputPixelType, InputImageDimension>;
  using ThresholdFilterType = BinaryThresholdImageFilter< InputImageType, TOutputImage >;
  using ErodeFilterType = BinaryErodeImageFilter< TOutputImage, TOutputImage, TKernel >;
  using DilateFilterType = BinaryDilateImageFilter< TOutputImage, TOutputImage, TKernel >;
  using SubtractFilterType = SubtractImageFilter< TOutputImage, TOutputImage, TOutputImage >;

protected:
  SelectNarrowbandImageFilter();
  ~SelectNarrowbandImageFilter() override = default;

  void GenerateData () override;

private:
  /* Member variables */
  SizeValueType m_Radius;
  typename ThresholdFilterType::Pointer m_ThresholdFilter;
  typename ErodeFilterType::Pointer m_ErodeFilter;
  typename DilateFilterType::Pointer m_DilateFilter;
  typename SubtractFilterType::Pointer m_SubtractFilter;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSelectNarrowbandImageFilter.hxx"
#endif

#endif // itkSelectNarrowbandImageFilter_h
