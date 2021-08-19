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
#ifndef itkSignedDitherImageFilter_h
#define itkSignedDitherImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{

/** \class SignedDitherImageFilter
 *
 * \brief Compute a signed dither of an image.
 *
 * The filter guaranttes the sign of the embedding does not change
 * when dithering.
 *
 * \author Bryce Besler
 * \ingroup Morphometry
 */
template <typename TImage>
class SignedDitherImageFilter : public ImageToImageFilter<TImage, TImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(SignedDitherImageFilter);

  static constexpr unsigned int ImageDimension = TImage::ImageDimension;

  /** Type defines */
  using PixelType = typename TImage::PixelType;
  using RealType = typename NumericTraits< PixelType >::RealType;
  using RegionType = typename TImage::RegionType;

  /** Standard class typedefs. */
  using Self = SignedDitherImageFilter<TImage>;
  using Superclass = ImageToImageFilter<TImage, TImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(SignedDitherImageFilter, ImageToImageFilter);

  /** Standard New macro. */
  itkNewMacro(Self);

  /** Set/Get Epsilon */
  itkSetMacro(Gamma, RealType);
  itkGetConstMacro(Gamma, RealType);

  /** Set if should divide gamma. */
  itkSetMacro(DivideGamma, bool);
  itkGetConstReferenceMacro(DivideGamma, bool);

  /** Set/Get fixed seed. */
  itkSetMacro(FixedSeed, bool);
  itkGetConstReferenceMacro(FixedSeed, bool);

protected:
  SignedDitherImageFilter();
  ~SignedDitherImageFilter() override = default;

  void PrintSelf(std::ostream & os, Indent indent) const override;

  void DynamicThreadedGenerateData(const RegionType & outputRegion) override;

private:
  /* Member variables */
  RealType m_Gamma;
  bool m_DivideGamma;
  bool m_FixedSeed;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSignedDitherImageFilter.hxx"
#endif

#endif // itkSignedDitherImageFilter_h
