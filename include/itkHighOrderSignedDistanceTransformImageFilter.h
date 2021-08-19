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
#ifndef itkHighOrderSignedDistanceTransformImageFilter_h
#define itkHighOrderSignedDistanceTransformImageFilter_h

#include "itkImageToImageFilter.h"

#include <queue>

namespace itk
{

/** \class HighOrderSignedDistanceTransformImageFilter
 *
 * \brief Signed distance transform of a binary signal with increased smoothness.
 *
 * 
 *
 * \author Bryce Besler
 * \ingroup Morphometry
 */
template <typename TInputImage, typename TOutputImage>
class HighOrderSignedDistanceTransformImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(HighOrderSignedDistanceTransformImageFilter);

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
  using Self = HighOrderSignedDistanceTransformImageFilter<InputImageType, TOutputImage>;
  using Superclass = ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(HighOrderSignedDistanceTransformImageFilter, ImageToImageFilter);

  /** Standard New macro. */
  itkNewMacro(Self);

  /** Heap type */
  using IndexType = typename OutputImageType::IndexType;
  using HeapElement = std::pair< RealType, IndexType >;
  using PriorityQueueType = std::priority_queue<HeapElement, std::vector<HeapElement>, std::greater< HeapElement >>;

protected:
  HighOrderSignedDistanceTransformImageFilter();
  ~HighOrderSignedDistanceTransformImageFilter() override = default;

  // Override since the filter needs all the data for the algorithm.
  void GenerateInputRequestedRegion() override;

  // Override since the filter produces the entire dataset.
  void EnlargeOutputRequestedRegion(DataObject * output) override;

  void GenerateData () override;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkHighOrderSignedDistanceTransformImageFilter.hxx"
#endif

#endif // itkHighOrderSignedDistanceTransformImageFilter_h
