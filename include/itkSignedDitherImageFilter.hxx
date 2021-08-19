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
#ifndef itkSignedDitherImageFilter_hxx
#define itkSignedDitherImageFilter_hxx

#include "itkSignedDitherImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include <random>

namespace itk
{

template <typename TImage>
SignedDitherImageFilter<TImage>
::SignedDitherImageFilter() :
    Superclass(),
    m_DivideGamma(true),
    m_FixedSeed(false),
    m_Gamma(10.0f)
{}

template <typename TImage>
void
SignedDitherImageFilter<TImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Gamma: " << m_Gamma << std::endl;
  os << indent << "Divide Gamma: " << m_DivideGamma << std::endl;
  os << indent << "Fixed seed: " << m_FixedSeed << std::endl;
}

template <typename TImage>
void
SignedDitherImageFilter<TImage>
::DynamicThreadedGenerateData(const RegionType & outputRegion)
{
  /* Create iterators */
  TImage *      output = this->GetOutput();
  const TImage * input = this->GetInput();

  itk::ImageRegionConstIterator<TImage> in(input, outputRegion);
  itk::ImageRegionIterator<TImage>     out(output, outputRegion);

  /* Set gamma */
  RealType gamma = this->m_Gamma;
  if (this->m_DivideGamma){
    RealType averageSpacing = 0.;
    for (auto &s: this->GetInput()->GetSpacing()){
      averageSpacing += s;
    }
    averageSpacing /= this->GetInput()->GetSpacing().GetNumberOfComponents();
    gamma = averageSpacing / this->m_Gamma;
  }

  /* Setup random */
  RealType amplitude, noise;
  std::default_random_engine generator;
  if (this->m_FixedSeed) {
    generator.seed(0);
  } else {
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
  }
  std::uniform_real_distribution<RealType> distribution(-1.0, +1.0);

  /* Iterator */
  for (in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd() && !out.IsAtEnd(); ++in, ++out)
  {
    amplitude = std::min(gamma, std::abs(in.Get()));
    noise = distribution(generator);
    out.Set(in.Get() + amplitude*noise);
  }
}

} // end namespace itk

#endif // itkSignedDitherImageFilter_hxx
