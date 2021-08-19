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
#ifndef itkCleanLabelSurfaceImageFilter_h
#define itkCleanLabelSurfaceImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkImageRegionIterator.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include <queue>
#include <set>

namespace itk
{

/** \class CleanLabelSurfaceImageFilter
 *
 * \brief Clean a binary surface of shocks
 *
 * A shock is defined as a pixel where the left and right neighbours are of
 * the same value, but both different from the center pixel.
 *    I(x-h) = I(x+h) =/= I(x)
 * This corresponds to a one pixel thin line across the image.
 *
 * There is no obvious way to parallelize this filter, so it can take some time
 * to process. For larger images, it is recommended to 
 *
 * \author Bryce Besler
 * \ingroup Morphometry
 */
template <typename TInputImage>
class CleanLabelSurfaceImageFilter : public ImageToImageFilter<TInputImage, TInputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(CleanLabelSurfaceImageFilter);

  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;

  /** Type defines */
  using InputImageType = TInputImage;
  using InputPixelType = typename InputImageType::PixelType;

  /** Standard class typedefs. */
  using Self = CleanLabelSurfaceImageFilter<InputImageType>;
  using Superclass = ImageToImageFilter<InputImageType, InputImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(CleanLabelSurfaceImageFilter, ImageToImageFilter);

  /** Standard New macro. */
  itkNewMacro(Self);

  /** Queue types */
  using BoundaryConditionType = typename itk::ZeroFluxNeumannBoundaryCondition<InputImageType, InputImageType>;
  using NeighborhoodIteratorType = typename itk::ConstShapedNeighborhoodIterator<InputImageType, BoundaryConditionType>;
  using OffsetType = typename NeighborhoodIteratorType::OffsetType;
  using OffsetsType = std::vector< OffsetType >;
  using IndexType = typename InputImageType::IndexType;
  using QueueType = std::queue< IndexType, std::deque< IndexType > >;

  /** Label types */
  using SizeType = ::itk::SizeValueType;
  using LabelVectorType = std::set< InputPixelType >;

  /** NumberOfLabels macros */
  itkSetMacro(NumberOfLabels, SizeType);
  itkGetConstMacro(NumberOfLabels, SizeType);

  /** MaxNumberOfLabels macros */
  itkSetMacro(MaxNumberOfLabels, SizeType);
  itkGetConstMacro(MaxNumberOfLabels, SizeType);

  /** BackgroundLabel macros */
  itkSetMacro(BackgroundLabel, InputPixelType);
  itkGetConstMacro(BackgroundLabel, InputPixelType);

protected:
  SizeType       m_NumberOfLabels;
  SizeType       m_MaxNumberOfLabels;
  InputPixelType m_BackgroundLabel;

  CleanLabelSurfaceImageFilter();
  ~CleanLabelSurfaceImageFilter() override = default;

  /** Override since the filter needs all the data for the algorithm. */
  void GenerateInputRequestedRegion() override;

  /** Override since the filter produces the entire dataset. */
  void EnlargeOutputRequestedRegion(DataObject * output) override;

  void GenerateData () override;

  /** Helper functions */
  void ProcessLabel(InputPixelType label, TInputImage * image, const OffsetsType &offsets);
  void ProcessQueue(QueueType &in, QueueType &out, InputPixelType label, TInputImage * image, const OffsetsType &offsets);
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkCleanLabelSurfaceImageFilter.hxx"
#endif

#endif // itkCleanLabelSurfaceImageFilter_h
