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
#ifndef itkSelectNarrowbandImageFilter_hxx
#define itkSelectNarrowbandImageFilter_hxx

#include "itkSelectNarrowbandImageFilter.h"
#include "itkNumericTraits.h"

namespace itk
{

template <typename TInputImage, typename TOutputImage>
SelectNarrowbandImageFilter<TInputImage, TOutputImage>
::SelectNarrowbandImageFilter() :
  Superclass(),
  m_Radius(1)
{
  this->m_ThresholdFilter = ThresholdFilterType::New();
  this->m_ErodeFilter = ErodeFilterType::New();
  this->m_DilateFilter = DilateFilterType::New();
  this->m_SubtractFilter = SubtractFilterType::New();

  /* Threshold */
  this->m_ThresholdFilter->SetUpperThreshold(NumericTraits<InputPixelType>::Zero);
  this->m_ThresholdFilter->SetInsideValue(NumericTraits<OutputPixelType>::One);
  this->m_ThresholdFilter->SetOutsideValue(NumericTraits<OutputPixelType>::Zero);

  /* Dilate */
  this->m_DilateFilter->SetInput(this->m_ThresholdFilter->GetOutput());
  this->m_DilateFilter->SetForegroundValue(NumericTraits<OutputPixelType>::One);
  this->m_DilateFilter->SetBackgroundValue(NumericTraits<OutputPixelType>::Zero);

  /* Erode */
  this->m_ErodeFilter->SetInput(this->m_ThresholdFilter->GetOutput());
  this->m_ErodeFilter->SetForegroundValue(NumericTraits<OutputPixelType>::One);
  this->m_ErodeFilter->SetBackgroundValue(NumericTraits<OutputPixelType>::Zero);

  /* Subtract */
  this->m_SubtractFilter->SetInput1(this->m_DilateFilter->GetOutput());
  this->m_SubtractFilter->SetInput2(this->m_ErodeFilter->GetOutput());
}

template <typename TInputImage, typename TOutputImage>
void
SelectNarrowbandImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  /* Allocate the output */
  this->AllocateOutputs();

  /* Graft input */
  typename TInputImage::Pointer input = TInputImage::New();
  input->Graft( const_cast< TInputImage * >( this->GetInput() ));

  /* Input */
  this->m_ThresholdFilter->SetInput(input);

  /* Structuring eleent */
  TKernel structuringElement;
  structuringElement.SetRadius(m_Radius);
  structuringElement.CreateStructuringElement();

  this->m_DilateFilter->SetKernel(structuringElement);
  this->m_ErodeFilter->SetKernel(structuringElement);

  /* Graft output */
  this->m_SubtractFilter->GraftOutput( this->GetOutput() );
  this->m_SubtractFilter->Update();
  this->GraftOutput(this->m_SubtractFilter->GetOutput() );
}

} // end namespace itk

#endif // itkSelectNarrowbandImageFilter_hxx
