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
#ifndef itkHighOrderSignedFastSweeping_hxx
#define itkHighOrderSignedFastSweeping_hxx

#include "itkHighOrderSignedFastSweeping.h"

#include "itkImageRegionIterator.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkFiniteDifference.h"

#include "itkSignalFunctions.h"

#include <algorithm>

namespace itk
{

template <typename TInputImage>
HighOrderSignedFastSweeping<TInputImage>
::HighOrderSignedFastSweeping() :
    Superclass(),
    m_WENOEpsilon(1e-6),
    m_MaxError(1e-6),
    m_Error(0.),
    m_MaxIteration(10),
    m_Iteration(0),
    m_LowerShift(0.),
    m_UpperShift(0.),
    m_NarrowbandSize(NumericTraits<RealType>::max()),
    m_Epsilon(itk::Math::eps)
{}

template <typename TInputImage>
void
HighOrderSignedFastSweeping<TInputImage>
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
HighOrderSignedFastSweeping<TInputImage>
::EnlargeOutputRequestedRegion(DataObject * output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}

template <typename TInputImage>
void
HighOrderSignedFastSweeping<TInputImage>
::GenerateData()
{
  /* need to explicitly allocate outputs */
  Superclass::AllocateOutputs();

  /* Get Inputs */
  InputImageType *      output = this->GetOutput();
  const InputImageType * input = this->GetInput();

  /**
   * 1) Copy
   */
  itk::ImageRegionConstIterator<InputImageType> in(input, output->GetRequestedRegion());
  itk::ImageRegionIterator<InputImageType>     out(output, output->GetRequestedRegion());

  for (in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd() && !out.IsAtEnd(); ++in, ++out)
  {
    /* Copy through */
    out.Set(in.Get());
  }

  /**
   * 2) Signed Fast Sweeping Method
   */
  auto region = output->GetRequestedRegion();
  StencilsType stencils = this->GenerateStencils();
  this->m_Iteration = 0;
  this->m_Error = NumericTraits< RealType >::max();
  while ( (this->m_Error > this->m_MaxError) && (this->m_Iteration < this->m_MaxIteration) )
  {
    /* So very sorry. I couldn't take the time to implement
     * as fast sweeping iterator. Everything else works for
     * N-dimensions though ;)
     */

    /* Error */
    // this->m_LastError = this->m_ThisError;
    this->m_Error = 0;

    for (int x = 0; x < region.GetSize(0); ++x) {
      for (int y = 0; y < region.GetSize(1); ++y) {
        for (int z = 0; z < region.GetSize(2); ++z) {
          IndexType index = {{x, y, z}};
          RealType lastValue = output->GetPixel(index);
          RealType solution = this->SolveIndex(output, index, stencils);
          output->SetPixel(index, solution);
          this->m_Error += std::abs(lastValue - solution);
        }
      }
    }

    for (int x = 0; x < region.GetSize(0); ++x) {
      for (int y = 0; y < region.GetSize(1); ++y) {
        for (int z = region.GetSize(2)-1; z >= 0; --z) {
          IndexType index = {{x, y, z}};
          RealType lastValue = output->GetPixel(index);
          RealType solution = this->SolveIndex(output, index, stencils);
          output->SetPixel(index, solution);
          this->m_Error += std::abs(lastValue - solution);
        }
      }
    }

    for (int x = 0; x < region.GetSize(0); ++x) {
      for (int y = region.GetSize(1)-1; y >= 0; --y) {
        for (int z = 0; z < region.GetSize(2); ++z) {
          IndexType index = {{x, y, z}};
          RealType lastValue = output->GetPixel(index);
          RealType solution = this->SolveIndex(output, index, stencils);
          output->SetPixel(index, solution);
          this->m_Error += std::abs(lastValue - solution);
        }
      }
    }

    for (int x = 0; x < region.GetSize(0); ++x) {
      for (int y = region.GetSize(1)-1; y >= 0; --y) {
        for (int z = region.GetSize(2)-1; z >= 0; --z) {
          IndexType index = {{x, y, z}};
          RealType lastValue = output->GetPixel(index);
          RealType solution = this->SolveIndex(output, index, stencils);
          output->SetPixel(index, solution);
          this->m_Error += std::abs(lastValue - solution);
        }
      }
    }

    for (int x = region.GetSize(0)-1; x >= 0; --x) {
      for (int y = 0; y < region.GetSize(1); ++y) {
        for (int z = 0; z < region.GetSize(2); ++z) {
          IndexType index = {{x, y, z}};
          RealType lastValue = output->GetPixel(index);
          RealType solution = this->SolveIndex(output, index, stencils);
          output->SetPixel(index, solution);
          this->m_Error += std::abs(lastValue - solution);
        }
      }
    }

    for (int x = region.GetSize(0)-1; x >= 0; --x) {
      for (int y = 0; y < region.GetSize(1); ++y) {
        for (int z = region.GetSize(2)-1; z >= 0; --z) {
          IndexType index = {{x, y, z}};
          RealType lastValue = output->GetPixel(index);
          RealType solution = this->SolveIndex(output, index, stencils);
          output->SetPixel(index, solution);
          this->m_Error += std::abs(lastValue - solution);
        }
      }
    }

    for (int x = region.GetSize(0)-1; x >= 0; --x) {
      for (int y = region.GetSize(1)-1; y >= 0; --y) {
        for (int z = 0; z < region.GetSize(2); ++z) {
          IndexType index = {{x, y, z}};
          RealType lastValue = output->GetPixel(index);
          RealType solution = this->SolveIndex(output, index, stencils);
          output->SetPixel(index, solution);
          this->m_Error += std::abs(lastValue - solution);
        }
      }
    }

    this->m_LowerShift = +1*NumericTraits<RealType>::max();
    this->m_UpperShift = -1*NumericTraits<RealType>::max();
    for (int x = region.GetSize(0)-1; x >= 0; --x) {
      for (int y = region.GetSize(1)-1; y >= 0; --y) {
        for (int z = region.GetSize(2)-1; z >= 0; --z) {
          IndexType index = {{x, y, z}};
          RealType lastValue = output->GetPixel(index);
          RealType solution = this->SolveIndex(output, index, stencils);
          output->SetPixel(index, solution);
          this->m_Error += std::abs(lastValue - solution);

          if (solution > 0) {
            this->m_LowerShift = std::min(this->m_LowerShift, solution);
          } else {
            this->m_UpperShift = std::max(this->m_UpperShift, solution);
          }
        }
      }
    }

    /* Shift, avoid if not needed */
    RealType shift = 0.5*(this->m_UpperShift+this->m_LowerShift);
    if (std::abs(shift) > this->m_Epsilon) {
      for (int x = 0; x < region.GetSize(0); ++x) {
        for (int y = 0; y < region.GetSize(1); ++y) {
          for (int z = 0; z < region.GetSize(2); ++z) {
            IndexType index = {{x, y, z}};
            output->SetPixel(index, output->GetPixel(index) - shift);
          }
        }
      }
    }

    this->m_Error /= (double)(2*2*2*region.GetSize(0)*region.GetSize(1)*region.GetSize(2));
    this->m_Iteration++;
    // std::cout << "Iteration: " << this->m_Iteration << " Error: " << this->m_Error;
    // std::cout << " Upper: " << this->m_UpperShift << " Lower: " << this->m_LowerShift << " Shift: " << shift << std::endl;
  }
}

template <typename TInputImage>
typename HighOrderSignedFastSweeping<TInputImage>::RealType
HighOrderSignedFastSweeping<TInputImage>
::SolveIndex(typename InputImageType::Pointer image, IndexType &index, StencilsType &stencils)
{
  /* Outside Narrow band? */
  RealType lastValue = image->GetPixel(index);
  if (std::abs(lastValue) > this->m_NarrowbandSize){
    return lastValue;
  }
  RealType s = (lastValue > 0) ? +1. : -1.;

  /* Sample */
  VectorVectorPairs neighbours = this->SampleStencil(image, index, stencils);

  /* Solve */
  return this->SolveUpwindQuadratic(neighbours, s);
}

template <typename TInputImage>
typename HighOrderSignedFastSweeping<TInputImage>::VectorVectorPairs
HighOrderSignedFastSweeping<TInputImage>
::SampleStencil(typename InputImageType::Pointer image, IndexType &index, StencilsType &stencils)
{
  /* Constants and sample */
  auto spacing = image->GetSpacing();
  auto region = image->GetRequestedRegion();
  RealType lastValue = image->GetPixel(index);
  RealType s = (lastValue > 0) ? +1. : -1.;

  /* Compute sided differences */
  VectorVectorPairs neighbours;
  for (unsigned int i = 0; i < InputImageDimension; ++i)
  {
    /* Stencil: samples[0], samples[1], samples[2], samples[3], samples[4] */
    using SamplesType = itk::Vector< RealType, 2*Order+1 >;
    SamplesType samples;
    for (unsigned int j = 0; j < 2*Order+1; ++j) {
      auto thisIndex = index + stencils[i][j];
      if (region.IsInside(thisIndex)) {
        samples[j] = image->GetPixel(thisIndex);
      } else {
        // samples[j] = lastValue; /* World's cheapest NN extrapolater :) */
        samples[j] = lastValue + s * spacing[i] * std::abs((RealType)(j - Order)); /* World's cheapest linear extrapolater :) */
      }
    }

    /* Left: samples[0], samples[1], samples[2], samples[3] */
    RealType left = weno_negative_fifth_order<RealType>(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5], spacing[i], this->m_WENOEpsilon);
    left = lastValue - spacing[i]*left; /* Eqn 2.12 */

    /* Right: samples[1], samples[2], samples[3], samples[4] */
    RealType right = weno_positive_fifth_order<RealType>(samples[1], samples[2], samples[3], samples[4], samples[5], samples[6], spacing[i], this->m_WENOEpsilon);
    right = lastValue + spacing[i]*right; /* Eqn 2.12 */

    RealType neighbour = s*std::min(s*left, s*right);
    neighbours[i] = {neighbour, spacing[i]};
  }

  return neighbours;
}

template <typename TInputImage>
typename HighOrderSignedFastSweeping<TInputImage>::RealType
HighOrderSignedFastSweeping<TInputImage>
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

  for (unsigned int j = 0; j < InputImageDimension; ++j)
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
  solution = s*std::max(s*solution, this->m_Epsilon);

  return solution;
}

template <typename TInputImage>
typename HighOrderSignedFastSweeping<TInputImage>::StencilsType
HighOrderSignedFastSweeping<TInputImage>
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

} // end namespace itk

#endif // itkHighOrderSignedFastSweeping_hxx
