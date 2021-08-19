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
#ifndef itkHighOrderSignedDistanceTransformCorrectionImageFilter_h
#define itkHighOrderSignedDistanceTransformCorrectionImageFilter_h

#include "itkImageToImageFilter.h"

#include <queue>

namespace itk
{

/** \class HighOrderSignedDistanceTransformCorrectionImageFilter
 *
 * \brief Signed distance transform of a binary signal with increased smoothness.
 *
 * 
 *
 * \author Bryce Besler
 * \ingroup Morphometry
 */
template <typename TInputImage>
class HighOrderSignedDistanceTransformCorrectionImageFilter : public ImageToImageFilter<TInputImage, TInputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(HighOrderSignedDistanceTransformCorrectionImageFilter);

  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int Order = 3;

  /** Type defines */
  using InputImageType = TInputImage;
  using InputPixelType = typename InputImageType::PixelType;
  using RealType = typename NumericTraits< InputPixelType >::RealType;
  using OffsetType = typename InputImageType::OffsetType;

  /** Standard class typedefs. */
  using Self = HighOrderSignedDistanceTransformCorrectionImageFilter<InputImageType>;
  using Superclass = ImageToImageFilter<InputImageType, TInputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(HighOrderSignedDistanceTransformCorrectionImageFilter, ImageToImageFilter);

  /** Standard New macro. */
  itkNewMacro(Self);

  /** Heap type */
  using IndexType = typename TInputImage::IndexType;
  using HeapElement = std::pair< RealType, IndexType >;
  using PriorityQueueType = std::priority_queue<HeapElement, std::vector<HeapElement>, std::greater< HeapElement >>;
  // using PriorityQueueType = std::queue<HeapElement, std::vector<HeapElement>, std::greater< HeapElement >>;


  /** Visited type */
  using VisitedImageType = itk::Image< unsigned char, InputImageDimension >;
  using StencilType = std::vector< OffsetType >;

protected:
  HighOrderSignedDistanceTransformCorrectionImageFilter();
  ~HighOrderSignedDistanceTransformCorrectionImageFilter() override = default;

  // Override since the filter needs all the data for the algorithm.
  void GenerateInputRequestedRegion() override;

  // Override since the filter produces the entire dataset.
  void EnlargeOutputRequestedRegion(DataObject * output) override;

  void GenerateData () override;

  RealType SolveIndex(typename InputImageType::Pointer image, typename VisitedImageType::Pointer visistedImage, IndexType &index, std::vector< StencilType > &stencils);
  typename std::vector< StencilType > GenerateStencils();
  StencilType GetNeighbors();
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkHighOrderSignedDistanceTransformCorrectionImageFilter.hxx"
#endif

#endif // itkHighOrderSignedDistanceTransformCorrectionImageFilter_h
