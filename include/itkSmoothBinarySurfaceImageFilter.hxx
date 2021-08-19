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
#ifndef itkSmoothBinarySurfaceImageFilter_hxx
#define itkSmoothBinarySurfaceImageFilter_hxx

#include "itkSmoothBinarySurfaceImageFilter.h"
#include "itkImageScanlineIterator.h"
#include "itkShapedNeighborhoodIterator.h"

#include <queue>
#include <vector>

namespace itk
{

template <typename TInputImage>
SmoothBinarySurfaceImageFilter<TInputImage>
::SmoothBinarySurfaceImageFilter() :
    Superclass()
{}

template <typename TInputImage>
void
SmoothBinarySurfaceImageFilter<TInputImage>
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
SmoothBinarySurfaceImageFilter<TInputImage>
::EnlargeOutputRequestedRegion(DataObject * output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}

template <typename TInputImage>
void
SmoothBinarySurfaceImageFilter<TInputImage>
::DynamicThreadedGenerateData (const ImageRegionType &outputRegionForThread)
{
  /* Create iterators */
  TInputImage * output = this->GetOutput();
  const TInputImage * input = this->GetInput();

  ImageScanlineConstIterator<TInputImage> inputIt(input, outputRegionForThread);
  ImageScanlineIterator<TInputImage>      outputIt(output, outputRegionForThread);

  while (!inputIt.IsAtEnd())
  {
    while (!inputIt.IsAtEndOfLine())
    {
      outputIt.Set(inputIt.Get()/1.1);
      ++inputIt;
      ++outputIt;
    }

    inputIt.NextLine();
    outputIt.NextLine();
  }
}

template <typename TInputImage>
void
SmoothBinarySurfaceImageFilter<TInputImage>
::AfterThreadedGenerateData()
{
  TInputImage * output = this->GetOutput();
  const TInputImage * input = this->GetInput();

  typename TInputImage::Pointer temp_image = TInputImage::New();
  temp_image->SetRegions(output->GetRequestedRegion());
  temp_image->Allocate();
  temp_image->FillBuffer(itk::NumericTraits<typename TInputImage::PixelType>::Zero);

  /* Narrow band */
  RealType gamma = 3.0;
  // using QueueType = std::vector<IndexType>;
  std::vector<IndexType> current, next;

  ImageScanlineIterator<TInputImage> outputIt(output, output->GetRequestedRegion());
  ImageScanlineIterator<TInputImage> tempIt(temp_image, output->GetRequestedRegion());

  while (!outputIt.IsAtEnd())
  {
    while (!outputIt.IsAtEndOfLine())
    {
      if (std::abs(outputIt.Get()) < gamma) {
        next.push_back(outputIt.GetIndex());
      }
      tempIt.Set(outputIt.Get());
      ++outputIt;
      ++tempIt;
    }
    outputIt.NextLine();
    tempIt.NextLine();
  }

  std::cout << next.size() << std::endl;

  RealType size = 1.;
  for(auto &s: output->GetRequestedRegion().GetSize()){
    size *= s;
  }
  StencilsType stencils = this->GenerateStencils();
  auto spacing = output->GetSpacing();
  auto region = output->GetRequestedRegion();

  RealType error = NumericTraits<RealType>::max();
  unsigned long i = 0;
  while ( (error > 1e-12) && (i < 1000) ) {
    error = 0.;
    int l = 0;
    int t = 0;
    std::swap(current, next);
    while (!current.empty()) {
      IndexType index = current.back();
      current.pop_back();

      // RealType margin = std::tanh(input->GetPixel(index) / 3.0);
      RealType margin = input->GetPixel(index);
      RealType sign = (margin >= 0) ? (+1.) : (-1.);
      RealType lastSolution;
      if (i%2) {
        lastSolution = output->GetPixel(index);
      } else {
        lastSolution = temp_image->GetPixel(index);
      }

      RealType a = 0.;
      RealType b = 0.;
      RealType c = 0.;
      // −1/12 4/3 −5/2 4/3 −1/12
      // Vector< RealType, 2*Order > coefficients = {{-1./12., 4./3., 4./3., -1./12.}};
      // std::vector< RealType> coefficients = {{-1./12., 4./3., 4./3., -1./12.}};
      // 1/90 −3/20 3/2 −49/18 3/2 −3/20 1/90
      std::vector< RealType> coefficients = {{1./90., -3/20., 3./2., 3./2., -3./20., 1./90.}};
      for (unsigned int i = 0; i < InputImageDimension; ++i) {
        RealType q = 0;
        for (unsigned int j = 0; j < 2*Order; ++j) {
          IndexType thisIndex = index + stencils[i][j];
          RealType qq = region.IsInside(thisIndex) ? (output->GetPixel(thisIndex)) : (lastSolution);
          q += qq*coefficients[j];
        }

        // a += std::pow(-5./2./spacing[i]/spacing[i], 2);
        a += std::pow(-49./18./spacing[i]/spacing[i], 2);
        b += 2. * -5./2./spacing[i]/spacing[i] * q/spacing[i]/spacing[i];
        c += std::pow(q/spacing[i]/spacing[i], 2);
      }

      /* Discrim */
      RealType solution;
      RealType discrim = b*b - 4.*a*c;
      if (discrim < 0.0) {
        solution = -b / (2.*a);
        // std::cout << a << " " << b << " " << c << std::endl;
        // itkExceptionMacro(<< "Discriminator is negative...");
        t++;
      } else {
        RealType s1 = (-b + std::sqrt(discrim)) / (2.*a);
        RealType s2 = (-b - std::sqrt(discrim)) / (2.*a);
        solution = sign*std::max(sign*s1, sign*s2);
        l++;
      }

      /* Solution */
      // RealType s1 = (-b + std::sqrt(discrim)) / (2.*a);
      // RealType s2 = (-b - std::sqrt(discrim)) / (2.*a);
      // RealType solution = sign*std::max(sign*s1, sign*s2);
      // RealType solution = -b / (2.*a);
      RealType w = 0.95;
      RealType propSolution = solution;
      solution = (1.-w) * solution + w * lastSolution;
      
      if (sign > 0) {
        solution = std::max(solution, margin);
      } else {
        solution = std::min(solution, margin);
      }

      // output->SetPixel(index, solution);
      if (i%2) {
        temp_image->SetPixel(index, solution);
        // lastSolution = output->GetPixel(index);
      } else{
        output->SetPixel(index, solution);
        // lastSolution = temp_image->GetPixel(index);
      }

      // std::cout << lastSolution << " " << solution << " " << margin << std::endl;

      next.push_back(index);

      error += std::abs(lastSolution - propSolution);
      // error += std::abs(lastSolution - solution);
    }
    error /= (RealType)(size);
    i++;
    std::cout << error << std::endl;
    std::cout << l << " " << t << std::endl;
  }
}

template <typename TInputImage>
typename SmoothBinarySurfaceImageFilter<TInputImage>::StencilsType
SmoothBinarySurfaceImageFilter<TInputImage>
::GenerateStencils()
{
  StencilsType stencils;
  for (unsigned int i = 0; i < InputImageDimension; ++i)
  {
    StencilType stencil;
    unsigned int k = 0;
    for (int j = -1*(int)Order; j <= (int)Order; ++j)
    {
      if (j == 0) {
        k = 1;
        continue;
      }
      OffsetType offset;
      offset.Fill(0);
      offset.SetElement(i, j);
      stencil[j + Order - k] = offset;
    }
    stencils[i] = stencil;
  }
  return stencils;
}

} // end namespace itk

#endif // itkSmoothBinarySurfaceImageFilter_hxx
