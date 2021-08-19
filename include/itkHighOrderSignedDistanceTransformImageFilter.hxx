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
#ifndef itkHighOrderSignedDistanceTransformImageFilter_hxx
#define itkHighOrderSignedDistanceTransformImageFilter_hxx

#include "itkHighOrderSignedDistanceTransformImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkFiniteDifference.h"

#include <algorithm>

namespace itk
{

template <typename TInputImage, typename TOutputImage>
HighOrderSignedDistanceTransformImageFilter<TInputImage, TOutputImage>
::HighOrderSignedDistanceTransformImageFilter() :
    Superclass()
{}

template <typename TInputImage, typename TOutputImage>
void
HighOrderSignedDistanceTransformImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
  if (this->GetInput())
  {
    auto * image = const_cast<InputImageType *>(this->GetInput());
    image->SetRequestedRegionToLargestPossibleRegion();
  }
}

template <typename TInputImage, typename TOutputImage>
void
HighOrderSignedDistanceTransformImageFilter<TInputImage, TOutputImage>
::EnlargeOutputRequestedRegion(DataObject * output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}

template <typename TInputImage, typename TOutputImage>
void
HighOrderSignedDistanceTransformImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  /* need to explicitly allocate outputs */
  Superclass::AllocateOutputs();

  /* Get Inputs */
  OutputImageType *      output = this->GetOutput();
  const InputImageType * input = this->GetInput();

  /* Setup heaps */
  PriorityQueueType store;
  PriorityQueueType process;

  /**
   * 1) Initialize
   */
  using BoundaryConditionType = typename itk::ZeroFluxNeumannBoundaryCondition<InputImageType, OutputImageType>;
  using NeighborhoodType = typename itk::ConstNeighborhoodIterator<InputImageType, BoundaryConditionType>;
  typename NeighborhoodType::RadiusType radius;
  radius.Fill(Order);

  NeighborhoodType in(radius, input, output->GetRequestedRegion());
  itk::ImageRegionIterator<OutputImageType> out(output, output->GetRequestedRegion());

  for (in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd() && !out.IsAtEnd(); ++in, ++out)
  {
    /* Center pixel and index */
    typename InputImageType::PointType x;
    input->TransformIndexToPhysicalPoint(in.GetIndex(), x);
    auto center_value = in.GetCenterPixel();

    /* Check if in narrow band */
    OutputPixelType distance = NumericTraits< OutputPixelType >::max();
    bool inNarrowBand = false;
    for (unsigned int i = 0; i < in.Size(); ++i)
    {
      if (in.GetPixel(i) != center_value)
      {
        inNarrowBand = true;
        typename InputImageType::PointType y;
        input->TransformIndexToPhysicalPoint(in.GetIndex(i), y);
        distance = std::min(distance, x.EuclideanDistanceTo(y));
      }
    }

    /* Set output */
    if (inNarrowBand)
    {
      if (center_value == 1)
      {
        distance = -1*distance/2.;
      }
      else
      {
        distance = distance/2.;
      }
      out.Set(distance);
      store.push({std::abs(distance), in.GetIndex()});
    }
    else
    {
      if (center_value == 1)
      {
        // out.Set(-1*NumericTraits< OutputPixelType >::max());
        out.Set(-100);
      }
      else
      {
        // out.Set(NumericTraits< OutputPixelType >::max());
        out.Set(100);
      }
    }
  }

  /**
   * 2) Solve Narrow Band
   */

  /* Setup neighbourhood */
  using StencilType = std::vector< typename NeighborhoodType::OffsetType >;
  std::vector< StencilType > stencils;
  for (unsigned int i = 0; i < OutputImageDimension; ++i)
  {
    StencilType stencil;
    for (int j = -1*(int)Order; j <= (int)Order; ++j)
    {
      typename NeighborhoodType::OffsetType offset;
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
  unsigned int max_iterations = 20;
  while (maxError > 0.1 and iterations < max_iterations)
  {
    std::swap(store, process);
    maxError = 0.;
    iterations++;

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
          RealType value = (region.IsInside(thisIndex)) ? output->GetPixel(thisIndex) : s*NumericTraits< RealType >::max();
          samples.push_back(value);
        }

        /* WENO sided differences */
        RealType left = weno_negative_fifth_order<RealType>(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5], spacing[i], epsilon);
        left = lastValue - spacing[i]*left; /* Eqn 2.12 */
        RealType right = weno_positive_fifth_order<RealType>(samples[1], samples[2], samples[3], samples[4], samples[5], samples[6], spacing[i], epsilon);
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

      /* Sort in absolute increasing order */
      std::sort(neighbours.begin(), neighbours.end());
      if (lastValue <= 0) {std::reverse(neighbours.begin(), neighbours.end());}

      /* Initial solution is very large */
      solution = s*NumericTraits< RealType >::max();

      // for (unsigned int j = 0; j < 1; ++j)
      for (unsigned int j = 1; j < neighbours.size(); ++j)
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
        solution = (-b + s*std::sqrt(discrim)) / (2.*a);
      }

      /* Compute */
      solution = s * std::max(s*solution, epsilon);

      /* Update */
      output->SetPixel(index, solution);
      store.push({std::abs(solution), index});
      maxError = std::max(maxError, std::abs(lastValue - solution));
    }

    std::cout << iterations << " " << maxError << std::endl;
  }

  std::cout << iterations << std::endl;


  /**
   * 3) Signed Fast Marching Method
   */
}

} // end namespace itk

#endif // itkHighOrderSignedDistanceTransformImageFilter_hxx
