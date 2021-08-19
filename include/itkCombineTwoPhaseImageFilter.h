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
#ifndef itkCombineTwoPhaseImageFilter_h
#define itkCombineTwoPhaseImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkNumericTraits.h"

#include "itkHeavisideImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"

namespace itk
{

/** \class CombineTwoPhaseImageFilter
 *
 * \brief Create a two phase image
 * 
 *  J = rho_1 [1 - H(-phi)] + rho_2 [H(-phi)]
 *
 * \author Bryce Besler
 * \ingroup Morphometry
 */
template <typename TImage>
class CombineTwoPhaseImageFilter : public ImageToImageFilter<TImage, TImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(CombineTwoPhaseImageFilter);

  static constexpr unsigned int ImageDimension = TImage::ImageDimension;

  /** Type defines */
  using PixelType = typename TImage::PixelType;
  using RealType = typename NumericTraits<PixelType>::RealType;

  /** Standard class typedefs. */
  using Self = CombineTwoPhaseImageFilter<TImage>;
  using Superclass = ImageToImageFilter<TImage, TImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(CombineTwoPhaseImageFilter, ImageToImageFilter);

  /** Standard New macro. */
  itkNewMacro(Self);

  /** Approximations to the Heaviside */
  itkSetMacro(Approximation, SmoothApproximationType);
  itkGetConstMacro(Approximation, SmoothApproximationType);

  void SetApproximationToTanh() {
    this->SetApproximation(Tanh);
  }

  void SetApproximationToSin() {
    this->SetApproximation(Sin);
  }

  /** Static method to convert enum to value */
  static std::string GetApproximationAsString(SmoothApproximationType approximation);

  /** Set/Get Epsilon */
  itkSetMacro(Epsilon, RealType);
  itkGetConstMacro(Epsilon, RealType);

  /** Set/Get Rho */
  itkSetMacro(Rho1, RealType);
  itkGetConstMacro(Rho1, RealType);
  itkSetMacro(Rho2, RealType);
  itkGetConstMacro(Rho2, RealType);

  /** Filter types */
  using HeavisideFilterType = HeavisideImageFilter< TImage >;
  using SubtractFilterType = SubtractImageFilter< TImage, TImage, TImage >;
  using MultiplyFilterType = MultiplyImageFilter< TImage >;
  using AddFilterType = AddImageFilter< TImage >;

protected:
  CombineTwoPhaseImageFilter();
  ~CombineTwoPhaseImageFilter() override = default;

  void GenerateData () override;

private:
  /* Member variables */
  SmoothApproximationType m_Approximation;
  RealType m_Epsilon;
  RealType m_Rho1;
  RealType m_Rho2;

  /* Filters */
  typename MultiplyFilterType::Pointer m_InverseFilter;
  typename HeavisideFilterType::Pointer m_HeavisideFilter;
  typename SubtractFilterType::Pointer m_SubtractFilter;
  typename MultiplyFilterType::Pointer m_MultiplyFilter1;
  typename MultiplyFilterType::Pointer m_MultiplyFilter2;
  typename AddFilterType::Pointer m_AddFilter;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkCombineTwoPhaseImageFilter.hxx"
#endif

#endif // itkCombineTwoPhaseImageFilter_h
