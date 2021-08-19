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
#ifndef itkDiracDeltaImageFilter_h
#define itkDiracDeltaImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkSignalFunctions.h"

namespace itk
{

/** \class DiracDeltaImageFilter
 *
 * \brief Compute the DiracDelta smooth approximation of an image
 *
 * The smooth approximations available are:
 *  Tanh:
 *      y = 1 / (epsilon * (cosh(2x/gamma) + 1))
 *  Sin:
 *      y = 0,                                         |x| > epsilon
 *          1/(2 * gamma) * (1 + cos(pi x /gamma)),    otherwise
 *
 * The primary difference is that the `tanh` approximation has an infinte
 * domain while the `sin` approximation has a finite domain.
 *
 * It is recommended to use floating point types with this filter.
 *
 * \author Bryce Besler
 * \ingroup Morphometry
 */
template <typename TInputImage, typename TOutputImage = TInputImage>
class DiracDeltaImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(DiracDeltaImageFilter);

  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

  /** Type defines */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using InputPixelType = typename InputImageType::PixelType;
  using OutputPixelType = typename OutputImageType::PixelType;
  using RealType = typename NumericTraits< InputPixelType >::RealType;

  /** Standard class typedefs. */
  using Self = DiracDeltaImageFilter<InputImageType, OutputImageType>;
  using Superclass = ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(DiracDeltaImageFilter, ImageToImageFilter);

  /** Standard New macro. */
  itkNewMacro(Self);

  /** Approximations to the DiracDelta */
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

protected:
  DiracDeltaImageFilter();
  ~DiracDeltaImageFilter() override = default;

  void PrintSelf(std::ostream & os, Indent indent) const override;

  using OutputRegionType = typename OutputImageType::RegionType;

  void DynamicThreadedGenerateData(const OutputRegionType & outputRegion) override;

private:
  /* Member variables */
  SmoothApproximationType m_Approximation;
  RealType m_Epsilon;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkDiracDeltaImageFilter.hxx"
#endif

#endif // itkDiracDeltaImageFilter_h
