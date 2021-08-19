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
#ifndef itkHeavisideImageFilter_h
#define itkHeavisideImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkSignalFunctions.h"

namespace itk
{

/** \class HeavisideImageFilter
 *
 * \brief Compute the Heaviside smooth approximation of an image
 *
 * The smooth approximations available are:
 *  Tanh:
 *      y = 1/2 * (1 + tanh(x/epsilon))
 *  Sin:
 *              0,                                               x < -epsilon
 *      y =     1,                                               x > +epsilon
 *              1/2 * (1 + x/epsilon + 1/pi sin(pi x /epsilon)), otherwise
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
class HeavisideImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(HeavisideImageFilter);

  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

  /** Type defines */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using InputPixelType = typename InputImageType::PixelType;
  using OutputPixelType = typename OutputImageType::PixelType;
  using RealType = typename NumericTraits< InputPixelType >::RealType;

  /** Standard class typedefs. */
  using Self = HeavisideImageFilter<InputImageType, OutputImageType>;
  using Superclass = ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(HeavisideImageFilter, ImageToImageFilter);

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

protected:
  HeavisideImageFilter();
  ~HeavisideImageFilter() override = default;

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
#  include "itkHeavisideImageFilter.hxx"
#endif

#endif // itkHeavisideImageFilter_h
