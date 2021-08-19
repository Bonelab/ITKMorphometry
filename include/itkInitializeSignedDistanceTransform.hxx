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
#ifndef itkInitializeSignedDistanceTransform_hxx
#define itkInitializeSignedDistanceTransform_hxx

#include "itkInitializeSignedDistanceTransform.h"

namespace itk
{

template <typename TInputImage, typename TOutputImage>
InitializeSignedDistanceTransform<TInputImage, TOutputImage>
::InitializeSignedDistanceTransform() :
  Superclass(),
  m_Dither(false),
  m_Label(1)
{
  this->m_ThresholdFilter = ThresholdFilterType::New();
  this->m_SubtractFilter = SubtractFilterType::New();
  this->m_ForegroundSDTFilter = DistanceMapFilterType::New();
  this->m_BackgroundSDTFilter = DistanceMapFilterType::New();
  this->m_SubtractSDTFilter = SubtractSDTFilerType::New();
  this->m_MultiplyFilter = MultiplyFilterType::New();
  this->m_DitherFilter = DitherFilterType::New();

  /* Threshold */
  this->m_ThresholdFilter->SetInsideValue(1);
  this->m_ThresholdFilter->SetOutsideValue(0);

  /* Invert */
  this->m_SubtractFilter->SetInput2(this->m_ThresholdFilter->GetOutput());
  this->m_SubtractFilter->SetConstant1(1);

  /* SDT foreground */
  this->m_ForegroundSDTFilter->SetInput(this->m_ThresholdFilter->GetOutput());
  this->m_ForegroundSDTFilter->InsideIsPositiveOff();
  this->m_ForegroundSDTFilter->UseImageSpacingOn();
  this->m_ForegroundSDTFilter->SquaredDistanceOff();

  /* SDT background */
  this->m_BackgroundSDTFilter->SetInput(this->m_SubtractFilter->GetOutput());
  this->m_BackgroundSDTFilter->InsideIsPositiveOff();
  this->m_BackgroundSDTFilter->UseImageSpacingOn();
  this->m_BackgroundSDTFilter->SquaredDistanceOff();

  /* Add images */
  this->m_SubtractSDTFilter->SetInput1(this->m_ForegroundSDTFilter->GetOutput());
  this->m_SubtractSDTFilter->SetInput2(this->m_BackgroundSDTFilter->GetOutput());

  /* Multiply */
  this->m_MultiplyFilter->SetInput1(this->m_SubtractSDTFilter->GetOutput());
  this->m_MultiplyFilter->SetConstant2(0.5);
}

template <typename TInputImage, typename TOutputImage>
void
InitializeSignedDistanceTransform<TInputImage, TOutputImage>
::GenerateData()
{
  /* Allocate the output */
  this->AllocateOutputs();

  /* Graft input */
  typename TInputImage::Pointer input = TInputImage::New();
  input->Graft( const_cast< TInputImage * >( this->GetInput() ));

  /* Input */
  this->m_ThresholdFilter->SetInput(input);

  /* Constant */
  this->m_ThresholdFilter->SetLowerThreshold(this->GetLabel());
  this->m_ThresholdFilter->SetUpperThreshold(this->GetLabel());

  /* Graft w/ dither */
  if (this->m_Dither) {
    m_DitherFilter->SetInput(m_MultiplyFilter->GetOutput());

    /* Graft output */
    m_DitherFilter->GraftOutput( this->GetOutput() );
    m_DitherFilter->Update();
    this->GraftOutput( m_DitherFilter->GetOutput() );
  } else {
    /* Graft output */
    m_MultiplyFilter->GraftOutput( this->GetOutput() );
    m_MultiplyFilter->Update();
    this->GraftOutput( m_MultiplyFilter->GetOutput() );
  }
}

} // end namespace itk

#endif // itkInitializeSignedDistanceTransform_hxx
