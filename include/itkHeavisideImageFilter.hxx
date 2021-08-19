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
#ifndef itkHeavisideImageFilter_hxx
#define itkHeavisideImageFilter_hxx

#include "itkHeavisideImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{

template <typename TInputImage, typename TOutputImage>
HeavisideImageFilter<TInputImage, TOutputImage>
::HeavisideImageFilter() :
    Superclass(),
    m_Approximation(Tanh),
    m_Epsilon(1.0f)
{}


template <typename TInputImage, typename TOutputImage>
void
HeavisideImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Epsilon: " << m_Epsilon << std::endl;
  os << indent << "Approximation: " << this->GetApproximationAsString(this->m_Approximation) << std::endl;
}

template <typename TInputImage, typename TOutputImage>
std::string
HeavisideImageFilter<TInputImage, TOutputImage>
::GetApproximationAsString(SmoothApproximationType approximation)
{
  switch(approximation) {
    case SmoothApproximationType::Tanh: return "tanh"; break;
    case SmoothApproximationType::Sin: return "sin"; break;
    default: return "Bad approximation"; break;
  }
}

template <typename TInputImage, typename TOutputImage>
void
HeavisideImageFilter<TInputImage, TOutputImage>
::DynamicThreadedGenerateData(const OutputRegionType & outputRegion)
{
  /* Create iterators */
  OutputImageType *      output = this->GetOutput();
  const InputImageType * input = this->GetInput();

  itk::ImageRegionConstIterator<InputImageType> in(input, outputRegion);
  itk::ImageRegionIterator<OutputImageType>     out(output, outputRegion);

  /* Create a function pointer to the method */
  RealType (*approx_function)(RealType, RealType);
  switch(this->GetApproximation()) {
    case SmoothApproximationType::Tanh: approx_function = &heaviside_tanh<RealType>; break;
    case SmoothApproximationType::Sin: approx_function = &heaviside_sin<RealType>; break;
    default: itkExceptionMacro(<< "Unknown approximation " << this->GetApproximation());
  }

  /* Iterator */
  for (in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd() && !out.IsAtEnd(); ++in, ++out)
  {
    out.Set(approx_function(in.Get(), this->m_Epsilon));
  }
}

} // end namespace itk

#endif // itkHeavisideImageFilter_hxx
