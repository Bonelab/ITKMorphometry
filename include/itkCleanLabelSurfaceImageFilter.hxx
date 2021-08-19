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
#ifndef itkCleanLabelSurfaceImageFilter_hxx
#define itkCleanLabelSurfaceImageFilter_hxx

#include "itkCleanLabelSurfaceImageFilter.h"

namespace itk
{

template <typename TInputImage>
CleanLabelSurfaceImageFilter<TInputImage>
::CleanLabelSurfaceImageFilter() :
    Superclass(),
    m_NumberOfLabels(NumericTraits< InputPixelType >::Zero),
    m_BackgroundLabel(NumericTraits< InputPixelType >::Zero),
    m_MaxNumberOfLabels(256)
{}

template <typename TInputImage>
void
CleanLabelSurfaceImageFilter<TInputImage>
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
CleanLabelSurfaceImageFilter<TInputImage>
::EnlargeOutputRequestedRegion(DataObject * output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}

template <typename TInputImage>
void
CleanLabelSurfaceImageFilter<TInputImage>
::GenerateData()
{
  /* need to explicitly allocate outputs */
  Superclass::AllocateOutputs();

  /* Get Inputs */
  InputImageType *      output = this->GetOutput();
  const InputImageType * input = this->GetInput();

  /* Labels */
  LabelVectorType labels;
  this->m_NumberOfLabels = 0;

  /* Queues for surface elements. */
  QueueType foregroundIn, foregroundProcess, backgroundIn, backgroundProcess;

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

  /* Setup iterators */
  itk::ImageRegionIterator<InputImageType> in(input, output->GetRequestedRegion());
  itk::ImageRegionIterator<InputImageType> out(output, output->GetRequestedRegion());

  /**
   * 1) Copy through and find number of labels
   */
  for (in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd() && !out.IsAtEnd(); ++in, ++out)
  {
    /* Pass through */
    auto center_value = in.GetCenterPixel();
    out.Set(center_value);

    /* Check if already seen */
    if ( (center_value != this->m_BackgroundLabel) && (labels.count(center_value) == 0))
    {
      labels.insert(center_value);
      this->m_NumberOfLabels++;
    
      /* Check size */
      if (this->m_NumberOfLabels >= this->m_MaxNumberOfLabels){
        itkExceptionMacro(<< "Number of labels exceeds max " << this->m_MaxNumberOfLabels << ". Did you accidentally provide a grey scale image?");
      }
    }
  }

  /**
   * 2) Process each label
   */
  for (auto &label: labels)
  {
    this->ProcessLabel(label, output, offsets);
  }
}

template <typename TInputImage>
void
CleanLabelSurfaceImageFilter<TInputImage>
::ProcessLabel(InputPixelType label, TInputImage * image, const OffsetsType &offsets)
{
  /* Queues for surface elements. */
  QueueType foregroundIn, foregroundProcess, backgroundIn, backgroundProcess;

  /* Setup iterators */
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(1);
  NeighborhoodIteratorType it(radius, image, image->GetRequestedRegion());

  /**
   * 1) Find surface elements
   */
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    /* Check if shock */
    auto center_value = (image->GetCenterPixel() == label);
    for (unsigned int i = 0; i < InputImageDimension; ++i)
    {
      auto leftValue = (image->GetPixel(offsets[2*i]) == label);
      auto rightValue = (image->GetPixel(offsets[2*i+1]) == label);
      if ((leftValue == rightValue) && (leftValue != center_value))
      {
        if (center_value == 0){
          backgroundIn.push(image->GetIndex());
        } else {
          foregroundIn.push(image->GetIndex());
        }
        break;
      }
    }
  }

  /**
   * 2) Process queue
   */
  while (backgroundIn.size() + foregroundIn.size() > 0)
  {
    std::swap(foregroundIn, foregroundProcess);
    std::swap(backgroundIn, backgroundProcess);

    /* Foreground */
    this->ProcessQueue(foregroundProcess, foregroundIn, label, image, offsets);

    /* Background */
    this->ProcessQueue(backgroundProcess, backgroundIn, label, image, offsets);
  }
}

template <typename TInputImage>
void
CleanLabelSurfaceImageFilter<TInputImage>
::ProcessQueue(QueueType &in, QueueType &out, InputPixelType label, InputImageType * image, const OffsetsType &offsets)
{
  auto region = image->GetRequestedRegion();
  while (!in.empty())
  {
    /* Pop */
    IndexType index = in.front();
    in.pop();

    /* Check if shock */
    auto center = (image->GetPixel(index) == label);
    bool shock = false;
    for (unsigned int i = 0; i < InputImageDimension; i++)
    {
      /* Compute offsets */
      auto leftIndex = index + offsets[2*i];
      auto rightIndex = index + offsets[2*i+1];

      /* Get Pixels */
      auto left = (region.IsInside(leftIndex)) ? image->GetPixel(leftIndex) : center;
      left = (left == label);
      auto right = (region.IsInside(rightIndex)) ? image->GetPixel(rightIndex) : center;
      right = (right == label);

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
      /* Invert */
      if (center) {
        image->SetPixel(index, label);
      } else {
        image->SetPixel(index, this->m_BackgroundLabel);
      }

      /* Check neighbours */
      for (auto &offset: offsets)
      {
        auto newIndex = index + offset;
        if (region.IsInside(newIndex) && (center == image->GetPixel(newIndex)))
        {
          out.push(newIndex);
        }
      }
    }
  }
}

} // end namespace itk

#endif // itkCleanLabelSurfaceImageFilter_hxx
