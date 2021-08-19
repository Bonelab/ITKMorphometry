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
#ifndef itkInitializeTwoPhaseNarrowBandImageFilter_hxx
#define itkInitializeTwoPhaseNarrowBandImageFilter_hxx

#include "itkInitializeTwoPhaseNarrowBandImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNumericTraits.h"
#include "itkMath.h"

namespace itk
{

template <typename TGrayImage, typename TNarrowbandImage>
InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>
::InitializeTwoPhaseNarrowBandImageFilter() :
  m_SplineOrder(3),
  m_ConvergenceOrder(3),
  m_MaxIterations(100)
{
  /* Inputs */
  this->SetNumberOfRequiredInputs(2);

  /* Outputs */
  this->SetNumberOfRequiredOutputs(2);
  this->SetNthOutput(0, static_cast<TGrayImage *>(this->MakeOutput(0).GetPointer()));
  this->SetNthOutput(1, static_cast<TNarrowbandImage *>(this->MakeOutput(1).GetPointer()));

  /* Interpolator */
  m_Interpolator = InterpolateType::New();
}

template <typename TGrayImage, typename TNarrowbandImage>
void
InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>
::SetGrayInput(const TGrayImage * greyImage)
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(0, const_cast<TGrayImage *>(greyImage));
}

template <typename TGrayImage, typename TNarrowbandImage>
const TGrayImage *
InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>
::GetGrayInput()
{
  return static_cast<const TGrayImage *>(this->ProcessObject::GetInput(0));
}

template <typename TGrayImage, typename TNarrowbandImage>
void
InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>
::SetNarrowbandInput(const TNarrowbandImage * narrowbandImage)
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(1, const_cast<TNarrowbandImage *>(narrowbandImage));
}

template <typename TGrayImage, typename TNarrowbandImage>
const TNarrowbandImage *
InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>
::GetNarrowbandInput()
{
  return static_cast<const TNarrowbandImage *>(this->ProcessObject::GetInput(1));
}

template <typename TGrayImage, typename TNarrowbandImage>
typename InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>::DataObjectPointer
InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>
::MakeOutput(DataObjectPointerArraySizeType idx)
{
  if (idx == 1)
  {
    return TNarrowbandImage::New().GetPointer();
  }
  return Superclass::MakeOutput(idx);
}

template <typename TGrayImage, typename TNarrowbandImage>
TGrayImage *
InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>
::GetInitializedNarrowbandOutput()
{
  return dynamic_cast<TGrayImage *>(this->ProcessObject::GetOutput(0));
}

template <typename TGrayImage, typename TNarrowbandImage>
TNarrowbandImage *
InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>
::GetNarrowbandMaskOutput()
{
  return dynamic_cast<TNarrowbandImage *>(this->ProcessObject::GetOutput(1));
}

template <typename TGrayImage, typename TNarrowbandImage>
void
InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>
::BeforeThreadedGenerateData()
{
  /* Setup interpolator */
  itkDebugMacro("Setting up interpolator");
  m_Interpolator->SetSplineOrder(m_SplineOrder);
  m_Interpolator->SetInputImage(this->GetGrayInput());

  /* Setup Epsilon */
  m_Epsilon = NumericTraits<RealType>::max();
  for (auto &s: this->GetGrayInput()->GetSpacing()) {
    m_Epsilon = std::min(m_Epsilon, s);
  }
  itkDebugMacro("Computed Epsilon: " << m_Epsilon);
}

template <typename TGrayImage, typename TNarrowbandImage>
void
InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>
::AfterThreadedGenerateData()
{
  /* Disconnect interpolator */
  m_Interpolator->SetInputImage(nullptr);
}

template <typename TGrayImage, typename TNarrowbandImage>
void
InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>
::DynamicThreadedGenerateData(const RegionType & region)
{
  /* Create iterators */
  const TGrayImage * greyInput = this->GetGrayInput();
  const TNarrowbandImage * narrowbandInput = this->GetNarrowbandInput();
  TGrayImage * initializedOutput = this->GetInitializedNarrowbandOutput();
  TNarrowbandImage * narrowbandOutput = this->GetNarrowbandMaskOutput();

  itk::ImageRegionConstIterator<TGrayImage>       greyIn(greyInput, region);
  itk::ImageRegionConstIterator<TNarrowbandImage> nbIn(narrowbandInput, region);
  itk::ImageRegionIterator<TGrayImage>            initOut(initializedOutput, region);
  itk::ImageRegionIterator<TNarrowbandImage>      nbOut(narrowbandOutput, region);

  /* Iterator */
  PointType inputPoint, closestPoint;
  SizeValueType iterations;
  for ( greyIn.GoToBegin(), nbIn.GoToBegin(), initOut.GoToBegin(), nbOut.GoToBegin();
        !greyIn.IsAtEnd() && !nbIn.IsAtEnd() && !initOut.IsAtEnd() && !nbOut.IsAtEnd();
        ++greyIn, ++nbIn, ++initOut, ++nbOut
      )
  {
    /* Get sign */
    RealType sign = (greyIn.Get() > 0 ) ? (+1.) : (-1.);

    /* Outside NB? */
    if (nbIn.Get() == NumericTraits<NarrowbandPixelType>::Zero) {
      /* Set NB mask */
      nbOut.Set(NumericTraits<NarrowbandPixelType>::Zero);

      /* Set init value */
      initOut.Set(sign*10.*m_Epsilon);
      continue;
    }

    /* Inside NB */
    greyInput->TransformIndexToPhysicalPoint(greyIn.GetIndex(), inputPoint);
    this->FindClosestPoint(inputPoint, closestPoint, iterations);

    /* Set */
    if (iterations < m_MaxIterations) {
      nbOut.Set(NumericTraits<NarrowbandPixelType>::One);
      initOut.Set(sign*inputPoint.EuclideanDistanceTo(closestPoint));
    } else {
      // std::cout << "No converge" << std::endl << std::flush;
      /* Deselect NB voxel */
      nbOut.Set(NumericTraits<NarrowbandPixelType>::Zero);

      /* Assign placeholder */
      initOut.Set(sign*10.*m_Epsilon);
    }
  }
}

template <typename TGrayImage, typename TNarrowbandImage>
void
InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>
::FindClosestPoint(PointType &inputPoint, PointType &closestPoint, SizeValueType& iterations)
{
  /* Initialize */
  closestPoint = inputPoint;
  iterations = 0;

  /* Run */
  // this->FindClosestPointFirstOrder(closestPoint, iterations);
  this->FindClosestPointPerpendicular(closestPoint, iterations);
}

// template <typename TGrayImage, typename TNarrowbandImage>
// void
// InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>
// ::FindClosestPointFirstOrder(PointType &point, SizeValueType& iterations)
// {
//   /* Initialize */
//   GrayPixelType value = 1.;
//   GrayPixelType norm = 0.;
//   RealType stepSize = m_Epsilon;
//   iterations = 0;

//   while ((std::abs(value) > std::pow(m_Epsilon, m_ConvergenceOrder) * norm) && (iterations < m_MaxIterations))
//   {
//     /* RK4 */
//     VectorType K1 = VelocityUpdate(point, value, norm);
//     VectorType K2 = VelocityUpdate(point + stepSize/2.*K1, value, norm);
//     VectorType K3 = VelocityUpdate(point + stepSize/2.*K2, value, norm);
//     VectorType K4 = VelocityUpdate(point + stepSize*K3, value, norm);

//     /* Update */
//     point = point + stepSize/6. * (K1 + 2.*K2 + 2.*K3 + K4);
//     iterations++;
//   }
//   // std::cout << iterations << " " << value << " " << norm << std::endl << std::flush;
// }

template <typename TGrayImage, typename TNarrowbandImage>
void
InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>
::FindClosestPointFirstOrder(PointType &point, SizeValueType& iterations)
{
  /* Initialize */
  GrayPixelType value = 1.;
  GrayPixelType norm = 0.;
  RealType stepSize = m_Epsilon;
  iterations = 0;

  while ((std::abs(value) > std::pow(m_Epsilon, m_ConvergenceOrder) * norm) && (iterations < m_MaxIterations))
  {
    /* RK4 */
    VectorType K1 = VelocityUpdate(point, value, norm);
    // VectorType K2 = VelocityUpdate(point + stepSize/2.*K1, value, norm);
    // VectorType K3 = VelocityUpdate(point + stepSize/2.*K2, value, norm);
    // VectorType K4 = VelocityUpdate(point + stepSize*K3, value, norm);

    /* Update */
    point = point + stepSize * (K1);
    iterations++;
  }
  // std::cout << "CP: " << iterations << " " << value << " " << norm << std::endl << std::flush;
}

template <typename TGrayImage, typename TNarrowbandImage>
typename InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>::VectorType
InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>
::VelocityUpdate(const PointType &point, GrayPixelType &value, GrayPixelType &norm)
{
  /* Value types */
  typename InterpolateType::CovariantVectorType grad;
  VectorType grad2;

  /* Interpolate */
  m_Interpolator->EvaluateValueAndDerivative(point, value, grad);
  norm = grad.Normalize();
  // norm = grad.GetNorm();

  /* Read into actual vector */
  for (unsigned int i = 0; i < GrayImageDimension; ++i)
  {
    grad2[i] = grad[i];
  }

  /* Update */
  return -1.*value / sqrt(value*value + norm*norm*m_Epsilon*m_Epsilon) * grad2;
}

template <typename TGrayImage, typename TNarrowbandImage>
void
InitializeTwoPhaseNarrowBandImageFilter<TGrayImage, TNarrowbandImage>
::FindClosestPointPerpendicular(PointType &point, SizeValueType& iterations)
{
  /* Initialize */
  PointType initialPoint = point;
  typename InterpolateType::CovariantVectorType grad;
  GrayPixelType value, gradNorm;
  GrayPixelType norm = 100. * std::pow(m_Epsilon, m_ConvergenceOrder);
  GrayPixelType d = 100. * std::pow(m_Epsilon, m_ConvergenceOrder);
  Vector< RealType, GrayImageDimension> normal, z, v;
  iterations = 0;
  SizeValueType cpIterations = 0;

  /* Initalize */
  this->FindClosestPointFirstOrder(point, cpIterations);

  /* Not converged, return early */
  if (cpIterations >= m_MaxIterations) {
    iterations = cpIterations;
    return;
  }

  while ((norm > std::pow(m_Epsilon, m_ConvergenceOrder)) && (iterations < m_MaxIterations))
  {
    /* Find normal */
    m_Interpolator->EvaluateValueAndDerivative(point, value, grad);
    gradNorm = grad.Normalize();
    for (unsigned int i = 0; i < GrayImageDimension; ++i)
    {
      normal[i] = grad[i];
    }

    /* Compute z */
    v = initialPoint.GetVectorFromOrigin() - point.GetVectorFromOrigin();
    z = v - (v * normal)*normal;

    /* Measure norm */
    norm = z.GetNorm();

    /* Update */
    // point = point + m_Epsilon*z;
    point = point + 0.5*z;
    
    /* Project */
    this->FindClosestPointFirstOrder(point, cpIterations);

    /* Not converged, return early */
    if (cpIterations >= m_MaxIterations) {
      iterations = cpIterations;
      return;
    }

    /* Update */
    iterations++;
    // std::cout << initialPoint << " " << iterations << " " << value << " " << norm << " " << std::pow(m_Epsilon, m_ConvergenceOrder) << std::endl << std::flush;
  }
  // std::cout << "PAR: " << iterations << " " << value << " " << norm << std::endl << std::flush;
}

} // end namespace itk

#endif // itkInitializeTwoPhaseNarrowBandImageFilter_hxx
