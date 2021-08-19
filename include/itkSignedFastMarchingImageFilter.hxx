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
#ifndef itkSignedFastMarchingImageFilter_hxx
#define itkSignedFastMarchingImageFilter_hxx

#include "itkSignedFastMarchingImageFilter.h"

#include "itkImageScanlineIterator.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkFiniteDifference.h"

#include <algorithm>

namespace itk
{

template <typename TInputImage>
SignedFastMarchingImageFilter<TInputImage>
::SignedFastMarchingImageFilter() :
    Superclass(),
    m_WENOEpsilon(1e-17),
    m_MaxError(-1),
    m_Error(0.),
    m_MaxIteration(10),
    m_Iteration(0),
    m_LowerShift(0.),
    m_UpperShift(0.),
    m_NarrowbandSize(NumericTraits<RealType>::max())
{}

template <typename TInputImage>
void
SignedFastMarchingImageFilter<TInputImage>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
  if (this->GetInput())
  {
    auto * image = const_cast<InputImageType *>(this->GetInput());
    image->SetRequestedRegionToLargestPossibleRegion();
  }
}

template <typename TInputImage>
void
SignedFastMarchingImageFilter<TInputImage>
::EnlargeOutputRequestedRegion(DataObject * output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}

template <typename TInputImage>
void
SignedFastMarchingImageFilter<TInputImage>
::GenerateData()
{
  /* need to explicitly allocate outputs */
  Superclass::AllocateOutputs();

  /* Get Inputs */
  InputImageType *      output = this->GetOutput();
  const InputImageType * input = this->GetInput();

  /* Setup queues */
  PriorityQueueType next, current;

  std::cout << "copy" << std::endl;

  /* Setup NB image */
  typename VisitedImageType::Pointer narrowbandImage = VisitedImageType::New();
  narrowbandImage->SetRegions(output->GetRequestedRegion());
  narrowbandImage->Allocate();
  narrowbandImage->FillBuffer(itk::NumericTraits<typename VisitedImageType::PixelType>::Zero);

  /**
   * 1) Copy
   */
  ImageScanlineConstIterator<TInputImage> inputIt(input, output->GetRequestedRegion());
  ImageScanlineIterator<TInputImage>      outputIt(output, output->GetRequestedRegion());
  ImageScanlineIterator<VisitedImageType> narrowbandIt(narrowbandImage, output->GetRequestedRegion());

  inputIt.GoToBegin();
  outputIt.GoToBegin();
  narrowbandIt.GoToBegin();
  while (!inputIt.IsAtEnd())
  {
    while (!inputIt.IsAtEndOfLine())
    {
      outputIt.Set(inputIt.Get());

      if (std::abs(inputIt.Get()) < this->m_NarrowbandSize) {
        next.push({std::abs(inputIt.Get()), inputIt.GetIndex()});
        narrowbandIt.Set(1);
      }

      ++inputIt;
      ++outputIt;
      ++narrowbandIt;
    }
    inputIt.NextLine();
    outputIt.NextLine();
    narrowbandIt.NextLine();
  }

  std::cout << " SFMM " << next.size() << std::endl;

  /**
   * 2) Signed Fast Marching Method
   */
  StencilsType stencils = this->GenerateStencils();
  StencilType neighbors = this->GetNeighbors();
  this->m_Iteration = 0;
  this->m_Error = NumericTraits< RealType >::max();
  while ( (this->m_Error > this->m_MaxError) && (this->m_Iteration < this->m_MaxIteration) )
  {
    /* SFMM */
    std::swap(next, current);
    this->SolveQueue(current, next, output, stencils, neighbors, narrowbandImage);

    /* Relevel */
    std::swap(next, current);
    this->Revel(current, next, output);

    this->m_Iteration++;
    std::cout << "Error: " << this->m_Error << " u: " << this->m_UpperShift << " l: " << this->m_LowerShift << " LogError: " << std::log(this->m_Error) << std::endl;
  }
}

template <typename TInputImage>
void
SignedFastMarchingImageFilter<TInputImage>
::SolveQueue(PriorityQueueType &current, PriorityQueueType &next, InputImageType * output, StencilsType &stencils, StencilType &neighbours, VisitedImageType * narrowbandImage)
{
  /* Setup */
  typename VisitedImageType::Pointer visistedImage = VisitedImageType::New();
  visistedImage->SetRegions(output->GetRequestedRegion());
  visistedImage->Allocate();
  visistedImage->FillBuffer(itk::NumericTraits<typename VisitedImageType::PixelType>::Zero);

  this->m_LowerShift = +1*NumericTraits<RealType>::max();
  this->m_UpperShift = -1*NumericTraits<RealType>::max();

  auto region = output->GetRequestedRegion();
  this->m_Error = 0;

  /* Solve */
  while (!current.empty()) {
    /* Pop */
    HeapElement element = current.top();
    current.pop();
    IndexType index = element.second;

    /* Visit */
    if (visistedImage->GetPixel(index) == 1){continue;}
    visistedImage->SetPixel(index, 1);
    next.push(element);

    RealType value = output->GetPixel(index);
    if (value > 0) {
      this->m_LowerShift = std::min(this->m_LowerShift, value);
    } else {
      this->m_UpperShift = std::max(this->m_UpperShift, value);
    }

    /* Visit neighbours */
    for (auto &neighbor: neighbours) {
      auto thisIndex = index + neighbor;
      if (region.IsInside(thisIndex) && (visistedImage->GetPixel(thisIndex) == 0) && (narrowbandImage->GetPixel(thisIndex) == 1)) {
        auto lastValue = output->GetPixel(thisIndex);
        auto solution = SampleAndSolve(output, thisIndex, stencils, narrowbandImage, visistedImage);
        output->SetPixel(thisIndex, solution);
        current.push({std::abs(solution), thisIndex});
        this->m_Error = std::max(this->m_Error, std::abs(lastValue - solution));
      }
    }
  }
}

template <typename TInputImage>
void
SignedFastMarchingImageFilter<TInputImage>
::Revel(PriorityQueueType &current, PriorityQueueType &next, InputImageType * output)
{
  /* Compute shift */
  RealType shift = 0.5*(this->m_LowerShift + this->m_UpperShift);
  if (std::abs(shift) < itk::Math::eps) {
    std::swap(current, next);
    return;
  }

  /* Solve */
  while (!current.empty()) {
    /* Pop */
    HeapElement element = current.top();
    current.pop();
    IndexType index = element.second;

    /* Visit */
    RealType value = output->GetPixel(index);
    output->SetPixel(index, value - shift);

    next.push({std::abs(value), index});
  }
}

template <typename TInputImage>
typename SignedFastMarchingImageFilter<TInputImage>::RealType
SignedFastMarchingImageFilter<TInputImage>
::SampleAndSolve(typename InputImageType::Pointer output, IndexType &index, StencilsType &stencils, VisitedImageType * narrowbandImage, VisitedImageType * visistedImage)
{
  /* 1) Sample */

  /* Constants and sample */
  auto spacing = output->GetSpacing();
  auto region = output->GetRequestedRegion();
  RealType lastValue = output->GetPixel(index);
  RealType s = (lastValue > 0) ? +1. : -1.;

  /* Compute sided differences */
  VectorVectorPairs neighbours;
  for (unsigned int i = 0; i < InputImageDimension; ++i)
  {
    /* Stencil: samples[0], samples[1], samples[2], samples[3], samples[4] */
    using SamplesType = itk::Vector< RealType, 2*Order+1 >;
    SamplesType samples;
    SamplesType valid;
    for (unsigned int j = 0; j < 2*Order+1; ++j) {
      auto thisIndex = index + stencils[i][j];
      if (region.IsInside(thisIndex) && (visistedImage->GetPixel(thisIndex) == 1) && (narrowbandImage->GetPixel(thisIndex) == 1)) {
        samples[j] = output->GetPixel(thisIndex);
        valid[j] = 1;
      } else {
        samples[j] = lastValue + s*spacing[i];
        valid[j] = 0;
      }
    }
    samples[Order] = lastValue;
    valid[Order] = 1;

    // /* Left: samples[0], samples[1], samples[2], samples[3] */
    // RealType left = s*NumericTraits< RealType >::max();
    // if (valid[1] == 1){
    //   if ( (valid[0] == 1) && (valid[3] == 1) ) {
    //     left = weno_negative_third_order<RealType>(samples[0], samples[1], samples[2], samples[3], spacing[i], this->m_WENOEpsilon);
    //     left = lastValue - spacing[i]*left; /* Eqn 2.12 */
    //   } else {
    //     left = samples[1];
    //   }
    // }

    // /* Right: samples[1], samples[2], samples[3], samples[4] */
    // RealType right = s*NumericTraits< RealType >::max();
    // if (valid[3] == 1){
    //   if ( (valid[1] == 1) && (valid[4] == 1) ) {
    //     right = weno_positive_third_order<RealType>(samples[1], samples[2], samples[3], samples[4], spacing[i], this->m_WENOEpsilon);
    //     right = lastValue + spacing[i]*right; /* Eqn 2.12 */
    //   } else {
    //     right = samples[3];
    //   }
    // }

    /* Left: samples[0], samples[1], samples[2], samples[3], samples[4], samples[5] */
    RealType left = s*NumericTraits< RealType >::max();
    if (valid[2] == 1){
      if ( (valid[1] == 1) && (valid[4] == 1) ) {
        if ( (valid[0] == 1) && (valid[5] == 1) ) {
          left = weno_negative_fifth_order<RealType>(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5], spacing[i], this->m_WENOEpsilon);
        } else {
          left = weno_negative_third_order<RealType>(samples[1], samples[2], samples[3], samples[4], spacing[i], this->m_WENOEpsilon);
        }
        left = lastValue - spacing[i]*left; /* Eqn 2.12 */
      } else {
        left = samples[2];
      }
    }

    /* Right: samples[1], samples[2], samples[3], samples[4], samples[5], samples[6] */
    RealType right = s*NumericTraits< RealType >::max();
    if (valid[4] == 1){
      if ( (valid[2] == 1) && (valid[5] == 1) ) {
        if ( (valid[0] == 1) && (valid[6] == 1) ) {
          right = weno_positive_fifth_order<RealType>(samples[1], samples[2], samples[3], samples[4], samples[5], samples[6], spacing[i], this->m_WENOEpsilon);
        } else {
          right = weno_positive_third_order<RealType>(samples[2], samples[3], samples[4], samples[5], spacing[i], this->m_WENOEpsilon);
        }
        right = lastValue + spacing[i]*right; /* Eqn 2.12 */
      } else {
        right = samples[4];
      }
    }

    if ( (left == right) && (left == s*NumericTraits< RealType >::max()) ) {continue;}

    RealType neighbour = s*std::min(s*left, s*right);
    // neighbours[i] = {neighbour, spacing[i]};
    neighbours.push_back({neighbour, spacing[i]});
  }

  /* 2) Solve */
  RealType solution = this->SolveUpwindQuadratic(neighbours, s);
  this->m_Error = std::abs(lastValue - solution);
  return solution;
}

template <typename TInputImage>
typename SignedFastMarchingImageFilter<TInputImage>::RealType
SignedFastMarchingImageFilter<TInputImage>
::SolveUpwindQuadratic(VectorVectorPairs &neighbours, RealType s)
{
  /* Solve quadratic */
  RealType a = 0.;
  RealType b = 0.;
  RealType c = -1.;
  RealType solution = 0.;
  RealType spaceFactor = 0.;
  RealType s1, s2;

  /* Sort in absolute increasing order */
  std::sort(neighbours.begin(), neighbours.end());
  if (s <= 0) {std::reverse(neighbours.begin(), neighbours.end());}

  /* Initial solution is very large */
  solution = s*NumericTraits< RealType >::max();

  for (unsigned int j = 0; j < neighbours.size(); ++j)
  {
    /* Break if between solutions */
    if (s*solution <= s*neighbours[j].first) {break;}

    /* Update this value */
    spaceFactor = 1. / (neighbours[j].second * neighbours[j].second);
    a += spaceFactor;
    b -= 2. * neighbours[j].first * spaceFactor;
    c += (neighbours[j].first * neighbours[j].first) * spaceFactor;

    /* Discrim */
    RealType discrim = b*b - 4.*a*c;
    if (discrim < 0.0)
    {
      std::cout << j << " " << a << " " << b << " " << c << " " << discrim << std::endl;
      itkExceptionMacro(<< "Discriminator is negative...");
    }

    /* Solution */
    s1 = (-b + std::sqrt(discrim)) / (2.*a);
    s2 = (-b - std::sqrt(discrim)) / (2.*a);
    solution = s*std::max(s*s1, s*s2);
  }

  /* Compute */
  solution = s * std::max(s*solution, itk::Math::eps);

  return solution;
}

template <typename TInputImage>
typename SignedFastMarchingImageFilter<TInputImage>::StencilsType
SignedFastMarchingImageFilter<TInputImage>
::GenerateStencils()
{
  StencilsType stencils;
  for (unsigned int i = 0; i < InputImageDimension; ++i)
  {
    StencilType stencil;
    for (int j = -1*(int)Order; j <= (int)Order; ++j)
    {
      OffsetType offset;
      offset.Fill(0);
      offset.SetElement(i, j);
      stencil[j + Order] = offset;
    }
    stencils[i] = stencil;
  }
  return stencils;
}

template <typename TInputImage>
typename SignedFastMarchingImageFilter<TInputImage>::StencilType
SignedFastMarchingImageFilter<TInputImage>
::GetNeighbors()
{
  StencilType neighbors;
  for (unsigned int i = 0; i < InputImageDimension; ++i)
  {
    /* Forward Element */
    OffsetType offset1;
    offset1.Fill(0);
    offset1.SetElement(i, +1);
    neighbors[2*i] = offset1;

    /* Backward Element */
    OffsetType offset2;
    offset2.Fill(0);
    offset2.SetElement(i, -1);
    neighbors[2*i+1] = offset2;
  }
  return neighbors;
}

} // end namespace itk

#endif // itkSignedFastMarchingImageFilter_hxx
