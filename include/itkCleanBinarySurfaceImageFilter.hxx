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
#ifndef itkCleanBinarySurfaceImageFilter_hxx
#define itkCleanBinarySurfaceImageFilter_hxx

#include "itkCleanBinarySurfaceImageFilter.h"
#include <cassert>

namespace itk
{

template <typename TInputImage>
CleanBinarySurfaceImageFilter<TInputImage>
::CleanBinarySurfaceImageFilter() :
    Superclass()
{}

template <typename TInputImage>
void
CleanBinarySurfaceImageFilter<TInputImage>
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
CleanBinarySurfaceImageFilter<TInputImage>
::EnlargeOutputRequestedRegion(DataObject * output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}


template <typename TInputImage>
void
CleanBinarySurfaceImageFilter<TInputImage>
::GenerateData()
{
  /* need to explicitly allocate outputs */
  Superclass::AllocateOutputs();

  /* Get Inputs */
  InputImageType *      output = this->GetOutput();
  const InputImageType * input = this->GetInput();

  /* Create blank image for checking */
  this->m_VisistedImage = VisitedImageType::New();
  this->m_VisistedImage->SetRegions(output->GetRequestedRegion());
  this->m_VisistedImage->Allocate();
  this->m_VisistedImage->FillBuffer(itk::NumericTraits<typename VisitedImageType::PixelType>::Zero);

  /* Queues for surface elements. */
  QueueType foregroundNext, foregroundProcess, backgroundNext, backgroundProcess;

  /* Setup iterators */
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(1);

  NeighborhoodIteratorType in(radius, input, output->GetRequestedRegion());
  itk::ImageRegionIterator<InputImageType>    out(output, output->GetRequestedRegion());

  /* Setup neighbourhood */
  OffsetsType offsets;
  for (unsigned int i = 0; i < InputImageDimension; ++i)
  {
    /* Forward Element */
    OffsetType offset1;
    offset1.Fill(0);
    offset1.SetElement(i, +1);
    offsets.push_back(offset1);

    /* Backward Element */
    OffsetType offset2;
    offset2.Fill(0);
    offset2.SetElement(i, -1);
    offsets.push_back(offset2);
  }

  /**
   * 1) Find surface elements
   */
  for (in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd() && !out.IsAtEnd(); ++in, ++out)
  {
    /* Pass through */
    auto center_value = 1 - static_cast<InputPixelType>(in.GetCenterPixel() == 0);
    out.Set(center_value);

    for (unsigned int i = 0; i < InputImageDimension; ++i)
    {
      auto leftValue = in.GetPixel(offsets[2*i]);
      auto rightValue = in.GetPixel(offsets[2*i+1]);
      if ((leftValue == rightValue) && (leftValue != center_value))
      {
        if (center_value == 0){
          backgroundNext.push(in.GetIndex());
        } else {
          foregroundNext.push(in.GetIndex());
        }
        break;
      }
    }
  }

  /**
   * 2) Process queue
   */
  while (backgroundNext.size() + foregroundNext.size() > 0)
  {
    std::swap(foregroundNext, foregroundProcess);
    std::swap(backgroundNext, backgroundProcess);

    this->ProcessQueue(backgroundProcess, foregroundNext, backgroundNext, output, offsets);
    this->ProcessQueue(foregroundProcess, foregroundNext, backgroundNext, output, offsets);
  }

  /* Clean up */
  this->m_VisistedImage = nullptr;
}

template <typename TInputImage>
void
CleanBinarySurfaceImageFilter<TInputImage>
::ProcessQueue(QueueType &in, QueueType &foregroundNext, QueueType &backgroundNext, InputImageType * image, const OffsetsType &offsets)
{
  auto region = image->GetRequestedRegion();
  while (!in.empty())
  {
    /* Pop */
    IndexType index = in.front();
    in.pop();

    /* Check if visisted */
    if (this->HaveVisisted(index)) continue;

    /* Check if shock */
    auto center = image->GetPixel(index);
    bool shock = false;
    for (unsigned int i = 0; i < InputImageDimension; i++)
    {
      /* Compute offsets */
      auto leftIndex = index + offsets[2*i];
      auto rightIndex = index + offsets[2*i+1];

      /* Get Pixels */
      auto left = (region.IsInside(leftIndex)) ? image->GetPixel(leftIndex) : center;
      auto right = (region.IsInside(rightIndex)) ? image->GetPixel(rightIndex) : center;

      /* Test for shock */
      if ((left == right) && (left != center))
      {
        shock = true;
        break;
      }
    }

    /* Process shock */
    if (shock)
    {
      /* Visit */
      this->MarkVisisted(index);

      /* Invert */
      image->SetPixel(index, 1 - center);

      /* Check neighbours */
      for (auto &offset: offsets)
      {
        auto newIndex = index + offset;
        if (region.IsInside(newIndex))
        {
          if (image->GetPixel(newIndex) == 0)
          {
            backgroundNext.push(newIndex);
          }
          else
          {
            foregroundNext.push(newIndex);
          }
        }
      }
    }
  }
}

template <typename TInputImage>
bool
CleanBinarySurfaceImageFilter<TInputImage>
::HaveVisisted(IndexType &index)
{
  if (this->m_VisistedImage->GetPixel(index) == 1) {
    return true;
  } else {
    return false;
  }
}

template <typename TInputImage>
void
CleanBinarySurfaceImageFilter<TInputImage>
::MarkVisisted(IndexType &index)
{
  this->m_VisistedImage->SetPixel(index, 1);
}

} // end namespace itk

#endif // itkCleanBinarySurfaceImageFilter_hxx
