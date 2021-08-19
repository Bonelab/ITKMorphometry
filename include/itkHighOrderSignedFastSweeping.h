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
#ifndef itkHighOrderSignedFastSweeping_h
#define itkHighOrderSignedFastSweeping_h

#include "itkImageToImageFilter.h"

namespace itk
{

/** \class HighOrderSignedFastSweeping
 *
 * \brief High order signed fast sweeping method.
 * 
 * Parameters:
 *  WENOEpsilon     WENO epsilon, typically 1e-6
 *  MaxError        Convergence criterion (average L1 norm across the 2^d sweeps)
 *  MaxIteration    Maximum number of iterations
 *  NarrowbandSize  Only update in a narrow band of specified size.
 * 
 * The input to this filter is a signed distance transform. If the signed distance
 * transform is computed from a sampled signal, it is highly recommended to slightly
 * dither the signal, or else convergence will be extremely slow due to a quantization
 * error in the signed distance transform of sampled signal.
 * 
 * We recommend at least two iterations when the signal is dithered. This produces very
 * reasonable results quickly. The use of narrowbanding is highly, highly recommended.
 * This converts the run time from scaling with the total number of voxels to only the
 * number of voxels on the surface, which is an incredible speed improvement for most
 * applications. Voxels outside the narrowband will be passed through. Since third
 * order WENO stencils are used, the narrow band should at least be selected at least as
 * large as three times the sample spacing.
 * 
 * The narrow band is computed as:
 *    abs(phi(x) ) < NarrowbandSize
 * So the narrow band will have thickness 2*NarrowbandSize. The parameter NarrowbandSize
 * can be interpreted as the distance from the surface over which to update.
 * 
 * Internal constants:
 *  Error           Error in current iteration
 *  Iteration       Current 
 *  LowerShift      Internal lower shift
 *  UpperShift      Internal upper shift
 * 
 * The embedding is shifted by subtracting 0.5*(LowerShift + UpperShift) at each iteration
 * to handle the non-uniqueness of the problem.
 *
 * This class currently only work for 3D images since I have not taken the time
 * to implement an N-dimensional fast sweeping iterator. However, everything else extends
 * N-dimensions.
 *
 * \author Bryce Besler
 * \ingroup Morphometry
 */
template <typename TInputImage>
class HighOrderSignedFastSweeping : public ImageToImageFilter<TInputImage, TInputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(HighOrderSignedFastSweeping);

  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int Order = 3;

  /** Type defines */
  using InputImageType = TInputImage;
  using InputPixelType = typename InputImageType::PixelType;
  using RealType = typename NumericTraits< InputPixelType >::RealType;
  using OffsetType = typename InputImageType::OffsetType;
  using IndexType = typename TInputImage::IndexType;
  using StencilType = itk::Vector< OffsetType, 2*Order +1 >;
  using StencilsType = itk::Vector< StencilType, InputImageDimension >;

  /** Standard class typedefs. */
  using Self = HighOrderSignedFastSweeping<InputImageType>;
  using Superclass = ImageToImageFilter<InputImageType, TInputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(HighOrderSignedFastSweeping, ImageToImageFilter);

  /** Standard New macro. */
  itkNewMacro(Self);

  /** Stencil pair macro.
   * These are fixed-size arrays which can be better optimized at compile time.
   */
  using VectorPairs = std::pair< RealType, RealType>;
  using VectorVectorPairs = itk::Vector< VectorPairs, InputImageDimension >;

  /** Set/Get WENOEpsilon */
  itkSetMacro(WENOEpsilon, RealType);
  itkGetConstMacro(WENOEpsilon, RealType);

  /** Set/Get MaxError */
  itkSetMacro(MaxError, RealType);
  itkGetConstMacro(MaxError, RealType);

  /** Get Error */
  itkGetConstMacro(Error, RealType);

  /** Set/Get MaxIteration */
  itkSetMacro(MaxIteration, unsigned int );
  itkGetConstMacro(MaxIteration, unsigned int);

  /** Get Iteration */
  itkGetConstMacro(Iteration, unsigned int);

  /** Get internal shift values */
  itkGetConstMacro(LowerShift, RealType);
  itkGetConstMacro(UpperShift, RealType);

  /** Set/Get narrow band value */
  itkSetMacro(NarrowbandSize, RealType);
  itkGetConstMacro(NarrowbandSize, RealType);

  /** Set/Get epsilon for preventing sign change */
  itkSetMacro(Epsilon, RealType);
  itkGetConstMacro(Epsilon, RealType);

protected:
  HighOrderSignedFastSweeping();
  ~HighOrderSignedFastSweeping() override = default;

  // Override since the filter needs all the data for the algorithm.
  void GenerateInputRequestedRegion() override;

  // Override since the filter produces the entire dataset.
  void EnlargeOutputRequestedRegion(DataObject * output) override;

  void GenerateData () override;

  // Helper classes for readability
  RealType SolveIndex(typename InputImageType::Pointer image, IndexType &index, StencilsType &stencils);
  StencilsType GenerateStencils();
  VectorVectorPairs SampleStencil(typename InputImageType::Pointer image, IndexType &index, StencilsType &stencils);
  RealType SolveUpwindQuadratic(VectorVectorPairs &neighbours, RealType s);

private:
  RealType m_WENOEpsilon;
  RealType m_MaxError;
  RealType m_Error;
  unsigned int m_MaxIteration;
  unsigned int m_Iteration;
  RealType m_LowerShift;
  RealType m_UpperShift;
  RealType m_NarrowbandSize;
  RealType m_Epsilon;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkHighOrderSignedFastSweeping.hxx"
#endif

#endif // itkHighOrderSignedFastSweeping_h
