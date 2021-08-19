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
#ifndef itkCombineTwoPhaseImageFilter_hxx
#define itkCombineTwoPhaseImageFilter_hxx

#include "itkCombineTwoPhaseImageFilter.h"
#include "itkNumericTraits.h"

namespace itk
{

template <typename TImage>
CombineTwoPhaseImageFilter<TImage>
::CombineTwoPhaseImageFilter() :
  Superclass(),
  m_Approximation(SmoothApproximationType::Sin),
  m_Epsilon(-1),
  m_Rho1(1),
  m_Rho2(0)
{
  this->m_InverseFilter = MultiplyFilterType::New();
  this->m_HeavisideFilter = HeavisideFilterType::New();
  this->m_SubtractFilter = SubtractFilterType::New();
  this->m_MultiplyFilter1 = MultiplyFilterType::New();
  this->m_MultiplyFilter2 = MultiplyFilterType::New();
  this->m_AddFilter = AddFilterType::New();

  /* Invert */
  this->m_InverseFilter->SetConstant1(-1 * NumericTraits<PixelType>::One);

  /* Heaviside */
  this->m_HeavisideFilter->SetInput(this->m_InverseFilter->GetOutput());

  /* Subtract */
  this->m_SubtractFilter->SetConstant1(NumericTraits<PixelType>::One);
  this->m_SubtractFilter->SetInput2(this->m_HeavisideFilter->GetOutput());

  /* Rho1 */
  this->m_MultiplyFilter1->SetInput1(this->m_SubtractFilter->GetOutput());

  /* Rho2 */
  this->m_MultiplyFilter2->SetInput1(this->m_HeavisideFilter->GetOutput());

  /* Combine */
  this->m_AddFilter->SetInput1(this->m_MultiplyFilter1->GetOutput());
  this->m_AddFilter->SetInput2(this->m_MultiplyFilter2->GetOutput());
}

template <typename TImage>
void
CombineTwoPhaseImageFilter<TImage>
::GenerateData()
{
  /* Allocate the output */
  this->AllocateOutputs();

  /* Graft input */
  typename TImage::Pointer input = TImage::New();
  input->Graft( const_cast< TImage * >( this->GetInput() ));

  /* Input */
  this->m_InverseFilter->SetInput2(input);

  /* Set parameters */
  this->m_MultiplyFilter1->SetConstant2(this->m_Rho1);
  this->m_MultiplyFilter2->SetConstant2(this->m_Rho2);

  this->m_HeavisideFilter->SetEpsilon(this->m_Epsilon);
  this->m_HeavisideFilter->SetApproximation(this->m_Approximation);

  /* Graft output */
  this->m_AddFilter->GraftOutput( this->GetOutput() );
  this->m_AddFilter->Update();
  this->GraftOutput(this->m_AddFilter->GetOutput() );
}

} // end namespace itk

#endif // itkCombineTwoPhaseImageFilter_hxx
