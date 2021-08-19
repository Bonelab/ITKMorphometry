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
#ifndef itkSignedFastMarchingImageFilter_h
#define itkSignedFastMarchingImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{

/** \class SignedFastMarchingImageFilter
 *
 * \brief High order signed fast marching method.
 *
 * \author Bryce Besler
 * \ingroup Morphometry
 */
template <typename TInputImage>
class SignedFastMarchingImageFilter : public ImageToImageFilter<TInputImage, TInputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(SignedFastMarchingImageFilter);

  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int Order = 3;

  /** Type defines */
  using InputImageType = TInputImage;
  using InputPixelType = typename InputImageType::PixelType;
  using RealType = typename NumericTraits< InputPixelType >::RealType;
  using OffsetType = typename InputImageType::OffsetType;
  using StencilType = itk::Vector< OffsetType, 2*Order +1 >;
  using StencilsType = itk::Vector< StencilType, InputImageDimension >;
  using RegionType = typename TInputImage::RegionType;

  /** Standard class typedefs. */
  using Self = SignedFastMarchingImageFilter<InputImageType>;
  using Superclass = ImageToImageFilter<InputImageType, TInputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(SignedFastMarchingImageFilter, ImageToImageFilter);

  /** Standard New macro. */
  itkNewMacro(Self);

  /** Stencil pair macro.
   * These are fixed-size arrays which can be better optimized at compile time.
   */
  using VectorPairs = std::pair< RealType, RealType>;
  // using VectorVectorPairs = itk::Vector< VectorPairs, InputImageDimension >;
  using VectorVectorPairs = std::vector< VectorPairs >;

  /** Heap type */
  using IndexType = typename TInputImage::IndexType;
  using HeapElement = std::pair< RealType, IndexType >;
  using PriorityQueueType = std::priority_queue<HeapElement, std::vector<HeapElement>, std::greater< HeapElement >>;

  /** Visited type */
  using VisitedImageType = itk::Image< unsigned char, InputImageDimension >;

  /** Set/Get WENOEpsilon */
  itkSetMacro(WENOEpsilon, RealType);
  itkGetConstMacro(WENOEpsilon, RealType);

  /** Set/Get MaxError */
  itkSetMacro(MaxError, RealType);
  itkGetConstMacro(MaxError, RealType);

  /** Get Error */
  itkGetConstMacro(Error, RealType);

  /** Set/Get MaxIteration */
  itkSetMacro(MaxIteration, unsigned int );
  itkGetConstMacro(MaxIteration, unsigned int);

  /** Get Iteration */
  itkGetConstMacro(Iteration, unsigned int);

  /** Get internal shift values */
  itkGetConstMacro(LowerShift, RealType);
  itkGetConstMacro(UpperShift, RealType);

  /** Set/Get narrow band value */
  itkSetMacro(NarrowbandSize, RealType);
  itkGetConstMacro(NarrowbandSize, RealType);

protected:
  SignedFastMarchingImageFilter();
  ~SignedFastMarchingImageFilter() override = default;

  // Override since the filter needs all the data for the algorithm.
  void GenerateInputRequestedRegion() override;
  void EnlargeOutputRequestedRegion(DataObject * output) override;

  void GenerateData () override;

  // Helper classes for readability
  void SolveQueue(PriorityQueueType &current, PriorityQueueType &next, InputImageType * output, StencilsType &stencils, StencilType &neighbours, VisitedImageType * narrowbandImage);
  RealType SampleAndSolve(typename InputImageType::Pointer output, IndexType &index, StencilsType &stencils, VisitedImageType * narrowbandImage,  VisitedImageType * visistedImage);
  RealType SolveUpwindQuadratic(VectorVectorPairs &neighbours, RealType s);
  void Revel(PriorityQueueType &current, PriorityQueueType &next, InputImageType * output);

  // Solver
  StencilsType GenerateStencils();
  StencilType GetNeighbors();

private:
  RealType m_WENOEpsilon;
  RealType m_MaxError;
  RealType m_Error;
  unsigned int m_MaxIteration;
  unsigned int m_Iteration;
  RealType m_LowerShift;
  RealType m_UpperShift;
  RealType m_NarrowbandSize;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSignedFastMarchingImageFilter.hxx"
#endif

#endif // itkSignedFastMarchingImageFilter_h
