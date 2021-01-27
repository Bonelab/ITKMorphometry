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
#ifndef itkCleanBinarySurfaceImageFilter_h
#define itkCleanBinarySurfaceImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkImageRegionIterator.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include <queue>

namespace itk
{

/** \class CleanBinarySurfaceImageFilter
 *
 * \brief Clean a binary surface of shocks
 *
 * A shock is defined as a pixel where the left and right neighbours are of
 * the same value, but both different from the center pixel.
 *    I(x-h) = I(x+h) =/= I(x)
 * This corresponds to a one pixel thin line across the image.
 *
 * This filter assumes the foreground label is '1' and the background label is '0'.
 * If the foreground label is not '1', the output will be '1'.
 *
 * \author Bryce Besler
 * \ingroup Morphometry
 */
template <typename TInputImage>
class CleanBinarySurfaceImageFilter : public ImageToImageFilter<TInputImage, TInputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(CleanBinarySurfaceImageFilter);

  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;

  /** Type defines */
  using InputImageType = TInputImage;
  using InputPixelType = typename InputImageType::PixelType;
  using VisitedImageType = itk::Image< unsigned char, InputImageDimension >;

  /** Standard class typedefs. */
  using Self = CleanBinarySurfaceImageFilter<InputImageType>;
  using Superclass = ImageToImageFilter<InputImageType, InputImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(CleanBinarySurfaceImageFilter, ImageToImageFilter);

  /** Standard New macro. */
  itkNewMacro(Self);

  /** Queue types */
  using BoundaryConditionType = typename itk::ZeroFluxNeumannBoundaryCondition<InputImageType, InputImageType>;
  using NeighborhoodIteratorType = typename itk::ConstShapedNeighborhoodIterator<InputImageType, BoundaryConditionType>;
  using OffsetType = typename NeighborhoodIteratorType::OffsetType;
  using OffsetsType = std::vector< OffsetType >;
  using IndexType = typename InputImageType::IndexType;
  using QueueType = std::queue< IndexType, std::deque< IndexType > >;

protected:
  CleanBinarySurfaceImageFilter();
  ~CleanBinarySurfaceImageFilter() override = default;

  /** Override since the filter needs all the data for the algorithm. */
  void GenerateInputRequestedRegion() override;

  /** Override since the filter produces the entire dataset. */
  void EnlargeOutputRequestedRegion(DataObject * output) override;

  void GenerateData () override;

  /** Helper functions */
  void ProcessQueue(QueueType &in, QueueType &foregroundNext, QueueType &backgroundNext, TInputImage * image, const OffsetsType &offsets);
  bool HaveVisisted(IndexType &index);
  void MarkVisisted(IndexType &index);

private:
 typename VisitedImageType::Pointer m_VisistedImage;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkCleanBinarySurfaceImageFilter.hxx"
#endif

#endif // itkCleanBinarySurfaceImageFilter_h
