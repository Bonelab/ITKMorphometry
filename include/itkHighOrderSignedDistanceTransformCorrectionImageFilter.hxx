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

// template <typename TInputImage>
// void
// HighOrderSignedDistanceTransformCorrectionImageFilter<TInputImage>
// ::GenerateData()
// {
//   /* need to explicitly allocate outputs */
//   Superclass::AllocateOutputs();

//   /* Get Inputs */
//   InputImageType *      output = this->GetOutput();
//   const InputImageType * input = this->GetInput();

//   /* Setup heaps */
//   PriorityQueueType queue, foregroundQueue, backgroundQueue;
//   IndexType temp;
//   temp.Fill(0);
//   HeapElement smallestForeground = {NumericTraits< RealType >::max(), temp};
//   HeapElement smallestBackground = {NumericTraits< RealType >::max(), temp};

//   /* Setup neighbourhood */
//   std::vector< StencilType > stencils = this->GenerateStencils();
//   StencilType neighbors = this->GetNeighbors();

//   /* Setup visited image */
//   typename VisitedImageType::Pointer visistedImage;

//   /* Setup narrowband image */
//   typename VisitedImageType::Pointer narrowbandImage = VisitedImageType::New();
//   narrowbandImage->SetRegions(output->GetRequestedRegion());
//   narrowbandImage->Allocate();
//   narrowbandImage->FillBuffer(itk::NumericTraits<typename VisitedImageType::PixelType>::Zero);

//   /* Setup narrow band value */
//   RealType narrowBandValue = 0.;
//   for (auto &s: output->GetSpacing()) {
//     narrowBandValue += s;
//   }
//   narrowBandValue *= 4. / output->GetSpacing().GetVectorDimension();
//   std::cout << narrowBandValue << std::endl;

//   /**
//    * 1) Copy
//    */
//   itk::ImageRegionConstIterator<InputImageType> in(input, output->GetRequestedRegion());
//   itk::ImageRegionIterator<InputImageType>     out(output, output->GetRequestedRegion());

//   for (in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd() && !out.IsAtEnd(); ++in, ++out)
//   {
//     /* Copy through */
//     out.Set(in.Get());

//     /* Get narrow band */
//     if (std::abs(in.Get()) < narrowBandValue) {
//       narrowbandImage->SetPixel(in.GetIndex(), 1);
//       if (in.Get() > 0) {
//         if (std::abs(in.Get()) < smallestBackground.first){
//           smallestBackground = {std::abs(in.Get()), in.GetIndex()};
//         }
//       } else {
//         if (std::abs(in.Get()) < smallestForeground.first){
//           smallestForeground = {std::abs(in.Get()), in.GetIndex()};
//         }
//       }
//     }
//   }

//   std::cout << "BKG: " << smallestBackground.first << " " << smallestBackground.second << std::endl;
//   std::cout << "FGN: " << smallestForeground.first << " " << smallestForeground.second << std::endl;

//   /**
//    * 2) Narrow Band
//    */
//   std::cout << "Narrow Band" << std::endl;
//   auto region = output->GetRequestedRegion();
//   unsigned int iterations = 0;
//   unsigned int maxIterations = 2;
//   RealType maxError = 1e-3;
//   RealType thisError = NumericTraits< RealType >::max();
//   while ( (thisError > maxError) && (iterations < maxIterations) )
//   {
//     /* Updates */
//     unsigned long long i = 0;
//     thisError = 0.;
//     iterations++;
//     visistedImage = VisitedImageType::New();
//     visistedImage->SetRegions(output->GetRequestedRegion());
//     visistedImage->Allocate();
//     visistedImage->FillBuffer(itk::NumericTraits<typename VisitedImageType::PixelType>::Zero);
//     RealType l = +1*NumericTraits<RealType>::max();
//     RealType u = -1*NumericTraits<RealType>::max();
//     foregroundQueue.push(smallestForeground);
//     backgroundQueue.push(smallestBackground);
//     smallestForeground = {NumericTraits< RealType >::max(), temp};
//     smallestBackground = {NumericTraits< RealType >::max(), temp};

//     /* Foreground */
//     while(!foregroundQueue.empty()) {
//       /* Pop */
//       HeapElement element = foregroundQueue.top();
//       foregroundQueue.pop();
//       IndexType index = element.second;

//       /* Check */
//       if (visistedImage->GetPixel(index) == 1) {continue;}

//       /* Mark visited */
//       visistedImage->SetPixel(index, 1);

//       /* Get value */
//       auto lastValue = output->GetPixel(index);

//       /* Set smallest */
//       if (std::abs(lastValue) < smallestForeground.first) {
//         smallestForeground = {std::abs(lastValue), index};
//       }

//       /* Store */
//       if (lastValue > 0) {
//         l = std::min(l, lastValue);
//       } else {
//         u = std::max(u, lastValue);
//       }

//       /* Solve rest */
//       for (auto &neighbor: neighbors) {
//         auto thisIndex = index + neighbor;
//         bool cond = region.IsInside(thisIndex) && (visistedImage->GetPixel(thisIndex) == 0) && (narrowbandImage->GetPixel(thisIndex) == 1) && (output->GetPixel(thisIndex) <= 0);
//         if (cond) {
//           auto lastValue = output->GetPixel(thisIndex);
//           auto solution = SolveIndex(output, narrowbandImage, thisIndex, stencils);
//           output->SetPixel(thisIndex, solution);
//           foregroundQueue.push({std::abs(solution), thisIndex});
//           // thisError = std::max(thisError, std::abs(lastValue - solution));
//           thisError += std::abs(lastValue - solution);
//           i+=1;
//         }
//       }
//     }

//     /* Background */
//     while(!backgroundQueue.empty()) {
//       /* Pop */
//       HeapElement element = backgroundQueue.top();
//       backgroundQueue.pop();
//       IndexType index = element.second;

//       /* Check */
//       if (visistedImage->GetPixel(index) == 1) {continue;}

//       /* Mark visited */
//       visistedImage->SetPixel(index, 1);

//       /* Get value */
//       auto lastValue = output->GetPixel(index);

//       /* Set smallest */
//       if (std::abs(lastValue) < smallestBackground.first) {
//         smallestBackground = {std::abs(lastValue), index};
//       }

//       /* Store */
//       if (lastValue > 0) {
//         l = std::min(l, lastValue);
//       } else {
//         u = std::max(u, lastValue);
//       }

//       /* Solve rest */
//       for (auto &neighbor: neighbors) {
//         auto thisIndex = index + neighbor;
//         bool cond = region.IsInside(thisIndex) && (visistedImage->GetPixel(thisIndex) == 0) && (narrowbandImage->GetPixel(thisIndex) == 1) && (output->GetPixel(thisIndex) > 0);
//         if (cond) {
//           auto lastValue = output->GetPixel(thisIndex);
//           auto solution = SolveIndex(output, narrowbandImage, thisIndex, stencils);
//           output->SetPixel(thisIndex, solution);
//           backgroundQueue.push({std::abs(solution), thisIndex});
//           // thisError = std::max(thisError, std::abs(lastValue - solution));
//           thisError += std::abs(lastValue - solution);
//           i+=1;
//         }
//       }
//     }

//     // /* March */
//     // while(!queue.empty())
//     // {
//     //   /* Pop */
//     //   HeapElement element = queue.top();
//     //   queue.pop();
//     //   IndexType index = element.second;

//     //   /* Check */
//     //   if (visistedImage->GetPixel(index) == 1) {continue;}

//     //   /* Mark visited */
//     //   visistedImage->SetPixel(index, 1);

//     //   /* Get value */
//     //   auto lastValue = output->GetPixel(index);

//     //   /* Set smallest */
//     //   if (std::abs(lastValue) < smallest.first) {
//     //     smallest = {std::abs(lastValue), index};
//     //   }

//     //   /* Store */
//     //   if (lastValue > 0) {
//     //     l = std::min(l, lastValue);
//     //   } else {
//     //     u = std::max(u, lastValue);
//     //   }

//     //   /* Solve rest */
//     //   for (auto &neighbor: neighbors) {
//     //     auto thisIndex = index + neighbor;
//     //     if (region.IsInside(thisIndex) && (visistedImage->GetPixel(thisIndex) == 0) && (narrowbandImage->GetPixel(thisIndex) == 1)) {
//     //       auto lastValue = output->GetPixel(thisIndex);
//     //       auto solution = SolveIndex(output, narrowbandImage, thisIndex, stencils);
//     //       output->SetPixel(thisIndex, solution);
//     //       queue.push({std::abs(solution), thisIndex});
//     //       // thisError = std::max(thisError, std::abs(lastValue - solution));
//     //       thisError += std::abs(lastValue - solution);
//     //       i+=1;
//     //     }
//     //   }
//     // }

//     if (i > 0) {
//       thisError /= (RealType)i;
//     }

//     // /* Relevel */
//     // if (std::abs(l - u) < 1e-6) {
//     //   visistedImage->FillBuffer(itk::NumericTraits<typename VisitedImageType::PixelType>::Zero);

//     //   while(!queue.empty())
//     //   {
//     //     /* Pop */
//     //     HeapElement element = queue.top();
//     //     queue.pop();
//     //     IndexType index = element.second;

//     //     /* Check */
//     //     if (visistedImage->GetPixel(index) == 1) {continue;}

//     //     /* Mark visited */
//     //     visistedImage->SetPixel(index, 1);

//     //     /* Relevel value */
//     //     auto lastValue = output->GetPixel(index) - 0.5*(l + u);
//     //     output->SetPixel(index, lastValue);

//     //     /* Set smallest */
//     //     if (std::abs(lastValue) < smallest.first) {
//     //       smallest = {std::abs(lastValue), in.GetIndex()};
//     //     }

//     //     /* Solve rest */
//     //     for (auto &neighbor: neighbors) {
//     //       auto thisIndex = index + neighbor;
//     //       if (region.IsInside(thisIndex) && (visistedImage->GetPixel(thisIndex) == 0) && (narrowbandImage->GetPixel(thisIndex) == 1)) {
//     //         queue.push({std::abs(output->GetPixel(thisIndex)), thisIndex});
//     //       }
//     //     }
//     //   }

//     //   queue.push(smallest);
//     //   smallest = {NumericTraits< RealType >::max(), temp};
//     // }
//     std::cout << "Iteration: " << iterations << " Error: " << thisError << " u: " << u << " l: " << l << std::endl;
//   }

//   /**
//    * 3) March
//    */
//   std::cout << "March" << std::endl;
//   while(!queue.empty())
//   {
//     /* Pop */
//     HeapElement element = queue.top();
//     queue.pop();
//     IndexType index = element.second;

//     /* Check */
//     if (visistedImage->GetPixel(index) == 1) {continue;}

//     /* Mark visited */
//     visistedImage->SetPixel(index, 1);

//     /* Get value */
//     auto lastValue = output->GetPixel(index);

//     /* Solve rest */
//     for (auto &neighbor: neighbors) {
//       auto thisIndex = index + neighbor;
//       if (region.IsInside(thisIndex) && (visistedImage->GetPixel(thisIndex) == 0)) {
//         auto lastValue = output->GetPixel(thisIndex);
//         auto solution = SolveIndex(output, visistedImage, thisIndex, stencils);
//         output->SetPixel(thisIndex, solution);
//         queue.push({std::abs(solution), thisIndex});
//       }
//     }
//   }
// }


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
  typename VisitedImageType::Pointer visistedImage;

  /* Setup narrowband image */
  typename VisitedImageType::Pointer narrowbandImage = VisitedImageType::New();
  narrowbandImage->SetRegions(output->GetRequestedRegion());
  narrowbandImage->Allocate();
  narrowbandImage->FillBuffer(itk::NumericTraits<typename VisitedImageType::PixelType>::Zero);

  /* Setup narrow band value */
  RealType narrowBandValue = 0.;
  for (auto &s: output->GetSpacing()) {
    narrowBandValue += s;
  }
  narrowBandValue *= 300000. / output->GetSpacing().GetVectorDimension();

  /**
   * 1) Copy
   */
  itk::ImageRegionConstIterator<InputImageType> in(input, output->GetRequestedRegion());
  itk::ImageRegionIterator<InputImageType>     out(output, output->GetRequestedRegion());

  for (in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd() && !out.IsAtEnd(); ++in, ++out)
  {
    /* Copy through */
    out.Set(in.Get());

    /* Get narrow band */
    if (std::abs(in.Get()) < narrowBandValue) {
      narrowbandImage->SetPixel(in.GetIndex(), 1);
      if (std::abs(in.Get()) < smallest.first) {
        smallest = {std::abs(in.Get()), in.GetIndex()};
      }
    }
  }

  /**
   * 2) Narrow Band
   */
  int i = 0;
  std::cout << "Narrow Band" << std::endl;
  auto region = output->GetRequestedRegion();
  unsigned int iterations = 0;
  unsigned int maxIterations = 100;
  RealType maxError = 1e-3;
  RealType thisError = NumericTraits< RealType >::max();
  while ( (thisError > maxError) && (iterations < maxIterations) )
  {
    /* Updates */
    unsigned long long i = 0;
    thisError = 0.;
    iterations++;
    visistedImage = VisitedImageType::New();
    visistedImage->SetRegions(output->GetRequestedRegion());
    visistedImage->Allocate();
    visistedImage->FillBuffer(itk::NumericTraits<typename VisitedImageType::PixelType>::Zero);
    RealType l = +1*NumericTraits<RealType>::max();
    RealType u = -1*NumericTraits<RealType>::max();
    queue.push(smallest);
    smallest = {NumericTraits< RealType >::max(), temp};

    /* March */
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
        smallest = {std::abs(lastValue), index};
      }

      /* Store */
      if (lastValue > 0) {
        l = std::min(l, lastValue);
      } else {
        u = std::max(u, lastValue);
      }

      /* Solve rest */
      for (auto &neighbor: neighbors) {
        auto thisIndex = index + neighbor;
        if (region.IsInside(thisIndex) && (visistedImage->GetPixel(thisIndex) == 0) && (narrowbandImage->GetPixel(thisIndex) == 1)) {
          auto lastValue = output->GetPixel(thisIndex);
          auto solution = SolveIndex(output, narrowbandImage, thisIndex, stencils);
          output->SetPixel(thisIndex, solution);
          queue.push({std::abs(solution), thisIndex});
          // thisError = std::max(thisError, std::abs(lastValue - solution));
          thisError += std::abs(lastValue - solution);
          i += 1;
        }
      }
    }

    if (i > 0) {
      thisError /= (RealType)i;
    }

    /* Relevel */
    if ( (std::abs(l - u) > 1e-9) && (l != +1*NumericTraits<RealType>::max()) && (u != -1*NumericTraits<RealType>::max()) ) {
      queue.push(smallest);
      smallest = {NumericTraits< RealType >::max(), temp};
      visistedImage = VisitedImageType::New();
      visistedImage->SetRegions(output->GetRequestedRegion());
      visistedImage->Allocate();
      visistedImage->FillBuffer(itk::NumericTraits<typename VisitedImageType::PixelType>::Zero);

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

        /* Relevel value */
        auto lastValue = output->GetPixel(index) - 0.5*(l + u);
        output->SetPixel(index, lastValue);

        /* Set smallest */
        if (std::abs(lastValue) < smallest.first) {
          smallest = {std::abs(lastValue), index};
        }

        /* Solve rest */
        for (auto &neighbor: neighbors) {
          auto thisIndex = index + neighbor;
          if (region.IsInside(thisIndex) && (visistedImage->GetPixel(thisIndex) == 0) && (narrowbandImage->GetPixel(thisIndex) == 1)) {
            queue.push({std::abs(output->GetPixel(thisIndex)), thisIndex});
          }
        }
      }
    }
    std::cout << "Iteration: " << iterations << " Error: " << thisError << " u: " << u << " l: " << l << std::endl;
  }

  /**
   * 3) March
   */
  unsigned int j = 0;
  queue.push(smallest);
  visistedImage = VisitedImageType::New();
  visistedImage->SetRegions(output->GetRequestedRegion());
  visistedImage->Allocate();
  visistedImage->FillBuffer(itk::NumericTraits<typename VisitedImageType::PixelType>::Zero);
  std::cout << "March " << queue.size() << std::endl;
  while(!queue.empty())
  {
    /* Pop */
    HeapElement element = queue.top();
    queue.pop();
    IndexType index = element.second;

    /* Check */
    if (visistedImage->GetPixel(index) == 1) {continue;}
    j++;

    /* Mark visited */
    visistedImage->SetPixel(index, 1);

    /* Get value */
    auto lastValue = output->GetPixel(index);

    /* Solve rest */
    for (auto &neighbor: neighbors) {
      auto thisIndex = index + neighbor;
      if (region.IsInside(thisIndex) && (visistedImage->GetPixel(thisIndex) == 0)) {
        auto lastValue = output->GetPixel(thisIndex);
        if (narrowbandImage->GetPixel(thisIndex) == 0) {
          auto solution = SolveIndex(output, visistedImage, thisIndex, stencils);
          output->SetPixel(thisIndex, solution);
          queue.push({std::abs(solution), thisIndex});
        }
        queue.push({std::abs(lastValue), thisIndex});
      }
    }
  }
  std::cout << j << std::endl;
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
  solution = s * std::max(s*solution, itk::Math::eps);

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
