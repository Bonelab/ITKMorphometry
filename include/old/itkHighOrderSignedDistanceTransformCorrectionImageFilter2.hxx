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
#ifndef itkHighOrderSignedDistanceTransformCorrectionImageFilter_hxx
#define itkHighOrderSignedDistanceTransformCorrectionImageFilter_hxx

#include "itkHighOrderSignedDistanceTransformCorrectionImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkFiniteDifference.h"

#include <algorithm>

namespace itk
{

template <typename TInputImage>
HighOrderSignedDistanceTransformCorrectionImageFilter<TInputImage>
::HighOrderSignedDistanceTransformCorrectionImageFilter() :
    Superclass()
{}

template <typename TInputImage>
void
HighOrderSignedDistanceTransformCorrectionImageFilter<TInputImage>
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
HighOrderSignedDistanceTransformCorrectionImageFilter<TInputImage>
::EnlargeOutputRequestedRegion(DataObject * output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}

template <typename TInputImage>
void
HighOrderSignedDistanceTransformCorrectionImageFilter<TInputImage>
::GenerateData()
{
  /* need to explicitly allocate outputs */
  Superclass::AllocateOutputs();

  /* Get Inputs */
  InputImageType *      output = this->GetOutput();
  const InputImageType * input = this->GetInput();

  /* Setup heaps */
  PriorityQueueType queue;
  IndexType temp;
  temp.Fill(0);
  HeapElement smallest = {NumericTraits< RealType >::max(), temp};

  /* Setup neighbourhood */
  std::vector< StencilType > stencils = this->GenerateStencils();
  StencilType neighbors = this->GetNeighbors();

  /* Setup visited image */
  typename VisitedImageType::Pointer visistedImage = VisitedImageType::New();
  visistedImage->SetRegions(output->GetRequestedRegion());
  visistedImage->Allocate();

  /**
   * 1) Copy
   */
  itk::ImageRegionConstIterator<InputImageType> in(input, output->GetRequestedRegion());
  itk::ImageRegionIterator<InputImageType>     out(output, output->GetRequestedRegion());

  for (in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd() && !out.IsAtEnd(); ++in, ++out)
  {
    /* Copy through */
    out.Set(in.Get());

    /* Get smallest */
    if (std::abs(in.Get()) < smallest.first) {
      smallest = {std::abs(in.Get()), in.GetIndex()};
    }
  }

  /**
   * 2) Narrow Band
   */


  /**
   * 3) March
   */

  auto region = output->GetRequestedRegion();
  unsigned int iterations = 0;
  unsigned int maxIterations = 1;
  RealType maxError = 0.01;
  RealType thisError = NumericTraits< RealType >::max();
  while ( (thisError > maxError) && (iterations < maxIterations) )
  {
    /* Updates */
    thisError = 0.;
    iterations++;
    visistedImage->FillBuffer(itk::NumericTraits<typename VisitedImageType::PixelType>::Zero);
    queue.push(smallest);

    /* March */
    std::cout << "March" << std::endl;
    while(!queue.empty())
    {
      /* Pop */
      HeapElement element = queue.top();
      queue.pop();
      IndexType index = element.second;

      /* Check */
      if (visistedImage->GetPixel(index) == 1) {continue;}

      /* Mark visited */
      visistedImage->SetPixel(index, 1);

      /* Get value */
      auto lastValue = output->GetPixel(index);

      /* Set smallest */
      if (std::abs(lastValue) < smallest.first) {
        smallest = {std::abs(lastValue), in.GetIndex()};
      }

      /* Solve rest */
      for (auto &neighbor: neighbors) {
        auto thisIndex = index + neighbor;
        if (region.IsInside(thisIndex) && (visistedImage->GetPixel(thisIndex) == 0)) {
          auto lastValue = output->GetPixel(thisIndex);
          auto solution = SolveIndex(output, visistedImage, thisIndex, stencils);
          output->SetPixel(thisIndex, solution);
          queue.push({std::abs(solution), thisIndex});
          thisError = std::max(thisError, std::abs(lastValue - solution));
        }
      }
    }
  }
}

template <typename TInputImage>
typename HighOrderSignedDistanceTransformCorrectionImageFilter<TInputImage>::RealType
HighOrderSignedDistanceTransformCorrectionImageFilter<TInputImage>
::SolveIndex(typename InputImageType::Pointer image, typename VisitedImageType::Pointer visistedImage, IndexType &index, std::vector< StencilType > &stencils)
{
  RealType epsilon = 1e-6;
  auto spacing = image->GetSpacing();
  auto region = image->GetRequestedRegion();

  /* Get */
  RealType lastValue = image->GetPixel(index);

  /* Sign */
  RealType s = (lastValue > 0) ? +1. : -1.;

  /* Compute sided differences */
  using VectorPairs = std::pair< RealType, RealType>;
  std::vector< VectorPairs > neighbours;
  for (unsigned int i = 0; i < stencils.size(); ++i)
  {
    /* Read neighbours */
    std::vector< RealType > samples;
    std::vector< unsigned char > valid;
    for (auto &offset: stencils[i])
    {
      auto thisIndex = index + offset;
      if (region.IsInside(thisIndex) && (visistedImage->GetPixel(thisIndex) == 1)) {
        samples.push_back(image->GetPixel(thisIndex));
        valid.push_back(1);
      } else {
        samples.push_back(0.);
        valid.push_back(0);
      }
    }
    samples[3] = lastValue;
    valid[3] = 1;

    /* Left: samples[0], samples[1], samples[2], samples[3], samples[4], samples[5] */
    RealType left = s*NumericTraits< RealType >::max();
    if (valid[2] == 1){
      if ( (valid[1] == 1) && (valid[4] == 1) ) {
        if ( (valid[0] == 1) && (valid[5] == 1) ) {
          left = weno_negative_fifth_order<RealType>(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5], spacing[i], epsilon);
        } else {
          left = weno_negative_third_order<RealType>(samples[1], samples[2], samples[3], samples[4], spacing[i], epsilon);
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
          right = weno_positive_fifth_order<RealType>(samples[1], samples[2], samples[3], samples[4], samples[5], samples[6], spacing[i], epsilon);
        } else {
          right = weno_positive_third_order<RealType>(samples[2], samples[3], samples[4], samples[5], spacing[i], epsilon);
        }
        right = lastValue + spacing[i]*right; /* Eqn 2.12 */
      } else {
        right = samples[4];
      }
    }

    if ( (left == right) && (left == s*NumericTraits< RealType >::max()) ) {return lastValue;}

    RealType neighbour = s*std::min(s*left, s*right);
    neighbours.push_back({neighbour, spacing[i]});
  }

  /* Solve quadratic */
  RealType a = 0.;
  RealType b = 0.;
  RealType c = -1.;
  RealType solution = 0.;
  RealType spaceFactor = 0.;
  RealType s1, s2;

  /* Sort in absolute increasing order */
  std::sort(neighbours.begin(), neighbours.end());
  if (lastValue <= 0) {std::reverse(neighbours.begin(), neighbours.end());}

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
  solution = s * std::max(s*solution, 1e-9);

  return solution;
}

template <typename TInputImage>
std::vector< typename HighOrderSignedDistanceTransformCorrectionImageFilter<TInputImage>::StencilType >
HighOrderSignedDistanceTransformCorrectionImageFilter<TInputImage>
::GenerateStencils()
{
  std::vector< StencilType > stencils;
  for (unsigned int i = 0; i < InputImageDimension; ++i)
  {
    StencilType stencil;
    for (int j = -1*(int)Order; j <= (int)Order; ++j)
    {
      OffsetType offset;
      offset.Fill(0);
      offset.SetElement(i, j);
      stencil.push_back(offset);
    }
    stencils.push_back(stencil);
  }
  return stencils;
}

template <typename TInputImage>
typename HighOrderSignedDistanceTransformCorrectionImageFilter<TInputImage>::StencilType
HighOrderSignedDistanceTransformCorrectionImageFilter<TInputImage>
::GetNeighbors()
{
  StencilType neighbors;
  for (unsigned int i = 0; i < InputImageDimension; ++i)
  {
    /* Forward Element */
    OffsetType offset1;
    offset1.Fill(0);
    offset1.SetElement(i, +1);
    neighbors.push_back(offset1);

    /* Backward Element */
    OffsetType offset2;
    offset2.Fill(0);
    offset2.SetElement(i, -1);
    neighbors.push_back(offset2);
  }
  return neighbors;
}

} // end namespace itk

#endif // itkHighOrderSignedDistanceTransformCorrectionImageFilter_hxx

// /* 2a) Solve smallest node */
// std::cout << "Smallest" << std::endl;
// IndexType index = smallest.second;
// auto lastValue = output->GetPixel(index);

// // Solve
// auto solution = SolveIndex(output, index, stencils);

// // Set
// output->SetPixel(index, solution);
// thisError = std::max(thisError, std::abs(lastValue - solution));
// smallest = {std::abs(solution), index};

// // Solve rest
// for (auto &neighbor: neighbors) {
//   auto thisIndex = index + neighbor;
//   if (region.IsInside(thisIndex) && (visistedImage->GetPixel(thisIndex) == 0)) {
//     auto lastValue = output->GetPixel(thisIndex);
//     auto solution = SolveIndex(output, thisIndex, stencils);
//     output->SetPixel(thisIndex, solution);
//     queue.push({std::abs(solution), thisIndex});
//     thisError = std::max(thisError, std::abs(lastValue - solution));
//   }
// }
