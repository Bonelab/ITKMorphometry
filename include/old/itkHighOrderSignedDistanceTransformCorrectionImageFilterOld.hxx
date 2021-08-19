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


/* Copy & NB */

/* Solve NB */

/* March */

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
  PriorityQueueType store;
  PriorityQueueType process;

  /**
   * 1) Copy
   */
  itk::ImageRegionConstIterator<InputImageType> in(input, output->GetRequestedRegion());
  itk::ImageRegionIterator<InputImageType>     out(output, output->GetRequestedRegion());


  RealType largeValue = 0.;
  for (in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd() && !out.IsAtEnd(); ++in, ++out)
  {
    /* Copy through */
    out.Set(in.Get());

    /* Get largest value */
    largeValue = std::max(largeValue, std::abs(in.Get()));

    /* Add to queue */
    store.push({std::abs(in.Get()), in.GetIndex()});
  }

  // You know, large!
  largeValue *= 2;

  std::cout << largeValue << std::endl;
  std::cout << store.size() << std::endl;

  /**
   * 2) Solve
   */

  /* Setup neighbourhood */
  using StencilType = std::vector< OffsetType >;
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

  RealType epsilon = 1e-6;
  auto spacing = output->GetSpacing();
  auto region = output->GetRequestedRegion();
  RealType maxError = NumericTraits< RealType >::max();
  unsigned int iterations = 0;
  unsigned int max_iterations = 1;
  while (maxError > 0.1 and iterations < max_iterations)
  {
    std::swap(store, process);
    maxError = 0.;
    iterations++;

    RealType l = +1*NumericTraits<RealType>::max();
    RealType u = -1*NumericTraits<RealType>::max();
    std::cout << l << " " << u << std::endl;

    while(!process.empty())
    {
      /* Pop */
      IndexType index = process.top().second;
      process.pop();

      /* Get */
      RealType lastValue = output->GetPixel(index);

      /* Sign */
      RealType s = (lastValue > 0) ? +1. : -1.;

      /* Compute sided differences */
      using VectorPairs = std::pair< RealType, RealType>;
      std::vector< VectorPairs > neighbours;
      for (unsigned int i = 0; i < stencils.size(); ++i)
      {
        /* Read neighbours */
        std::vector< RealType > samples;
        for (auto &offset: stencils[i])
        {
          auto thisIndex = index + offset;
          // RealType value = (region.IsInside(thisIndex)) ? output->GetPixel(thisIndex) : s*NumericTraits< RealType >::max();
          // RealType value = (region.IsInside(thisIndex)) ? output->GetPixel(thisIndex) : s*largeValue;
          RealType value = (region.IsInside(thisIndex)) ? output->GetPixel(thisIndex) : lastValue;
          samples.push_back(value);
        }

        /* WENO sided differences */
        // RealType left = weno_negative_fifth_order<RealType>(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5], spacing[i], epsilon);
        RealType left = weno_negative_third_order<RealType>(samples[1], samples[2], samples[3], samples[4], spacing[i], epsilon);
        left = lastValue - spacing[i]*left; /* Eqn 2.12 */
        // RealType right = weno_positive_fifth_order<RealType>(samples[1], samples[2], samples[3], samples[4], samples[5], samples[6], spacing[i], epsilon);
        RealType right = weno_positive_third_order<RealType>(samples[2], samples[3], samples[4], samples[5], spacing[i], epsilon);
        right = lastValue + spacing[i]*right; /* Eqn 2.12 */
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

      /* Store */
      if (s > 0) l = std::min(l, solution);
      if (s < 0) u = std::max(u, solution);

      /* Update */
      output->SetPixel(index, solution);
      store.push({std::abs(solution), index});
      maxError = std::max(maxError, std::abs(lastValue - solution));
    }

    RealType tt = 0.5*(l + u);

    std::cout << iterations << " " << maxError << std::endl;
    std::cout << u << " " << l << " " << tt << std::endl;

    std::swap(store, process);
    while(!process.empty())
    {
      /* Pop */
      IndexType index = process.top().second;
      process.pop();

      /* Shift */
      RealType lastValue = output->GetPixel(index);
      lastValue = lastValue - tt;

      /* Push */
      output->SetPixel(index, lastValue);
      store.push({std::abs(lastValue), index});
    }
  }

  std::cout << iterations << std::endl;
}

} // end namespace itk

#endif // itkHighOrderSignedDistanceTransformCorrectionImageFilter_hxx
