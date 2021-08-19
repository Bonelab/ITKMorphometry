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
#ifndef itkHighOrderSignedFastSweeping2_hxx
#define itkHighOrderSignedFastSweeping2_hxx

#include "itkHighOrderSignedFastSweeping2.h"

#include "itkImageRegionIterator.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkFiniteDifference.h"

#include "itkSignalFunctions.h"

#include <algorithm>

namespace itk
{

template <typename TInputImage, typename TBinaryImage>
HighOrderSignedFastSweeping2<TInputImage, TBinaryImage>
::HighOrderSignedFastSweeping2() :
    Superclass(),
    m_WENOEpsilon(1e-6),
    m_MaxError(1e-6),
    m_Error(0.),
    m_MaxIteration(10),
    m_Iteration(0),
    m_NarrowbandSize(NumericTraits<RealType>::max()),
    m_Epsilon(itk::Math::eps),
    m_Order(2)
{
  /* Inputs */
  this->SetNumberOfRequiredInputs(2);
}


template <typename TInputImage, typename TBinaryImage>
void
HighOrderSignedFastSweeping2<TInputImage, TBinaryImage>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
  if (this->GetInput())
  {
    auto * image = const_cast<InputImageType *>(this->GetInput());
    image->SetRequestedRegionToLargestPossibleRegion();
  }
}

template <typename TInputImage, typename TBinaryImage>
void
HighOrderSignedFastSweeping2<TInputImage, TBinaryImage>
::EnlargeOutputRequestedRegion(DataObject * output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}

template <typename TInputImage, typename TBinaryImage>
void
HighOrderSignedFastSweeping2<TInputImage, TBinaryImage>
::SetNarrowbandInput(const TInputImage * greyImage)
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(0, const_cast<TInputImage *>(greyImage));
}

template <typename TInputImage, typename TBinaryImage>
const TInputImage *
HighOrderSignedFastSweeping2<TInputImage, TBinaryImage>
::GetNarrowbandInput()
{
  return static_cast<const TInputImage *>(this->ProcessObject::GetInput(0));
}

template <typename TInputImage, typename TBinaryImage>
void
HighOrderSignedFastSweeping2<TInputImage, TBinaryImage>
::SetMaskInput(const TBinaryImage * narrowbandImage)
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(1, const_cast<TBinaryImage *>(narrowbandImage));
}

template <typename TInputImage, typename TBinaryImage>
const TBinaryImage *
HighOrderSignedFastSweeping2<TInputImage, TBinaryImage>
::GetMaskInput()
{
  return static_cast<const TBinaryImage *>(this->ProcessObject::GetInput(1));
}

template <typename TInputImage, typename TBinaryImage>
void
HighOrderSignedFastSweeping2<TInputImage, TBinaryImage>
::GenerateData()
{
  /* need to explicitly allocate outputs */
  Superclass::AllocateOutputs();

  /* Get Inputs */
  InputImageType *      output = this->GetOutput();
  const InputImageType * input = this->GetNarrowbandInput();
  const TBinaryImage *    mask = this->GetMaskInput();

  /**
   * 1) Copy
   */
  itk::ImageRegionConstIterator<InputImageType> in(input, output->GetRequestedRegion());
  itk::ImageRegionIterator<InputImageType>     out(output, output->GetRequestedRegion());

  for (in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd() && !out.IsAtEnd(); ++in, ++out)
  {
    /* Copy through */
    out.Set(std::abs(in.Get()));
  }

  /**
   * 2) Fast Sweeping Method
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
          RealType solution = this->SolveIndex(output, mask, index, stencils);
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
          RealType solution = this->SolveIndex(output, mask, index, stencils);
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
          RealType solution = this->SolveIndex(output, mask, index, stencils);
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
          RealType solution = this->SolveIndex(output, mask, index, stencils);
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
          RealType solution = this->SolveIndex(output, mask, index, stencils);
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
          RealType solution = this->SolveIndex(output, mask, index, stencils);
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
          RealType solution = this->SolveIndex(output, mask, index, stencils);
          output->SetPixel(index, solution);
          this->m_Error += std::abs(lastValue - solution);
        }
      }
    }

    for (int x = region.GetSize(0)-1; x >= 0; --x) {
      for (int y = region.GetSize(1)-1; y >= 0; --y) {
        for (int z = region.GetSize(2)-1; z >= 0; --z) {
          IndexType index = {{x, y, z}};
          RealType lastValue = output->GetPixel(index);
          RealType solution = this->SolveIndex(output, mask, index, stencils);
          output->SetPixel(index, solution);
          this->m_Error += std::abs(lastValue - solution);
        }
      }
    }

    this->m_Error /= (double)(2*2*2*region.GetSize(0)*region.GetSize(1)*region.GetSize(2));
    this->m_Iteration++;
    // std::cout << "Iteration: " << this->m_Iteration << " Error: " << this->m_Error << std::endl;
  }

  /**
   * 3) Sign
   */
  for (in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd() && !out.IsAtEnd(); ++in, ++out)
  {
    /* Get Sign */
    RealType s = in.Get() > 0 ? +1. : -1.;
    out.Set(s * out.Get());
  }
}

template <typename TInputImage, typename TBinaryImage>
typename HighOrderSignedFastSweeping2<TInputImage, TBinaryImage>::RealType
HighOrderSignedFastSweeping2<TInputImage, TBinaryImage>
::SolveIndex(typename InputImageType::Pointer image, typename TBinaryImage::ConstPointer mask, IndexType &index, StencilsType &stencils)
{
  /* Return if solved */
  if (mask->GetPixel(index) == 1) {
    return image->GetPixel(index);
  }

  /* Outside Narrow band? */
  RealType lastValue = image->GetPixel(index);
  if (std::abs(lastValue) > this->m_NarrowbandSize){
    return lastValue;
  }

  /* Sample */
  VectorVectorPairs neighbours = this->SampleStencil(image, index, stencils);

  /* Solve */
  return this->SolveUpwindQuadratic(neighbours);
}

template <typename TInputImage, typename TBinaryImage>
typename HighOrderSignedFastSweeping2<TInputImage, TBinaryImage>::VectorVectorPairs
HighOrderSignedFastSweeping2<TInputImage, TBinaryImage>
::SampleStencil(typename InputImageType::Pointer image, IndexType &index, StencilsType &stencils)
{
  /* Constants and sample */
  auto spacing = image->GetSpacing();
  auto region = image->GetRequestedRegion();
  RealType lastValue = image->GetPixel(index);

  /* Compute sided differences */
  VectorVectorPairs neighbours;
  for (unsigned int i = 0; i < InputImageDimension; ++i)
  {
    /* Stencil: samples[0], samples[1], samples[2], samples[3], samples[4] */
    using SamplesType = itk::Vector< RealType, 2*MaxOrder+1 >;
    SamplesType samples;
    for (unsigned int j = 0; j < 2*MaxOrder+1; ++j) {
      auto thisIndex = index + stencils[i][j];
      if (region.IsInside(thisIndex)) {
        samples[j] = image->GetPixel(thisIndex);
      } else {
        // samples[j] = lastValue; /* World's cheapest NN extrapolater :) */
        samples[j] = lastValue + spacing[i] * std::abs((RealType)(j - MaxOrder)); /* World's cheapest linear extrapolater :) */
      }
    }


    /* Left: samples[0], samples[1], samples[2], samples[3], samples[4], samples[5] */
    /* Right: samples[1], samples[2], samples[3], samples[4], samples[5], samples[6] */
    RealType left, right;
    switch(m_Order) {
      case 1:
        left = samples[2];
        right = samples[4];
        break;
      case 2:
        left = weno_negative_third_order<RealType>(samples[1], samples[2], samples[3], samples[4], spacing[i], this->m_WENOEpsilon);
        left = lastValue - spacing[i]*left; /* Eqn 2.12 */

        right = weno_positive_third_order<RealType>(samples[2], samples[3], samples[4], samples[5], spacing[i], this->m_WENOEpsilon);
        right = lastValue + spacing[i]*right; /* Eqn 2.12 */

        break;
      case 3:
        /* Left: samples[0], samples[1], samples[2], samples[3], samples[4], samples[5] */
        left = weno_negative_fifth_order<RealType>(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5], spacing[i], this->m_WENOEpsilon);
        left = lastValue - spacing[i]*left; /* Eqn 2.12 */

        /* Right: samples[1], samples[2], samples[3], samples[4], samples[5], samples[6] */
        right = weno_positive_fifth_order<RealType>(samples[1], samples[2], samples[3], samples[4], samples[5], samples[6], spacing[i], this->m_WENOEpsilon);
        right = lastValue + spacing[i]*right; /* Eqn 2.12 */

        break;
      default:
        itkExceptionMacro(<< "Unsupported order " << m_Order);
        break;
    }

    RealType neighbour = std::min(left, right);
    neighbours[i] = {neighbour, spacing[i]};
  }

  return neighbours;
}

template <typename TInputImage, typename TBinaryImage>
typename HighOrderSignedFastSweeping2<TInputImage, TBinaryImage>::RealType
HighOrderSignedFastSweeping2<TInputImage, TBinaryImage>
::SolveUpwindQuadratic(VectorVectorPairs &neighbours)
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

  /* Initial solution is very large */
  solution = NumericTraits< RealType >::max();

  for (unsigned int j = 0; j < InputImageDimension; ++j)
  {
    /* Break if between solutions */
    if (solution <= neighbours[j].first) {break;}

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
    solution = std::max(s1, s2);
  }

  return solution;
}

template <typename TInputImage, typename TBinaryImage>
typename HighOrderSignedFastSweeping2<TInputImage, TBinaryImage>::StencilsType
HighOrderSignedFastSweeping2<TInputImage, TBinaryImage>
::GenerateStencils()
{
  StencilsType stencils;
  for (unsigned int i = 0; i < InputImageDimension; ++i)
  {
    StencilType stencil;
    for (int j = -1*(int)MaxOrder; j <= (int)MaxOrder; ++j)
    {
      OffsetType offset;
      offset.Fill(0);
      offset.SetElement(i, j);
      stencil[j + MaxOrder] = offset;
    }
    stencils[i] = stencil;
  }
  return stencils;
}

} // end namespace itk

#endif // itkHighOrderSignedFastSweeping2_hxx
