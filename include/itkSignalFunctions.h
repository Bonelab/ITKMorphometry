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
#ifndef itkSignalFunctions_h
#define itkSignalFunctions_h

#include "itkMath.h"

namespace itk
{

typedef enum SmoothApproximation {
    Tanh = 1,
    Sin
} SmoothApproximationType;

/** Hyperbolic tan approximation to the Heaviside
 *  y = 1/2 * (1 + tanh(x/epsilon))
 */
template <typename T>
T heaviside_tanh(T x, T epsilon)
{
    return 1./2. * (1.0 + std::tanh(x/epsilon));
}

/** Hyperbolic tan approximation to the Dirac Delta
 *  y = 1 / (epsilon * (cosh(2x/gamma) + 1))
 */
template <typename T>
T dirac_delta_tanh(T x, T epsilon)
{
    T cosh = std::cosh(2.0 * x / epsilon);
    return 1.0 / (epsilon * (cosh + 1.0));
}

/** Sin approximation to the Heaviside
 *           0,                                              x < -epsilon
 *   y =     1,                                              x > +epsilon
 *           1/2 * (1 + x/epsilon + 1/pi sin(pi x /epsilon)),    otherwise
 */
template <typename T>
T heaviside_sin(T x, T epsilon)
{
    if (x < -epsilon)
    {
        return 0.;
    }
    else if (x > epsilon)
    {
        return 1.;
    }
    else
    {
        return 1./2. * (1. + x/epsilon + 1/Math::pi * std::sin(Math::pi * x /epsilon));
    }
}

/** Sin approximation to the Direc Delta
 *           0,                                          |x| > epsilon
 *   y =     1/(2 * gamma) * (1 + cos(pi x /gamma)),    otherwise
 */
template <typename T>
T dirac_delta_sin(T x, T epsilon)
{
    if (Math::abs(x) > epsilon)
    {
        return 0.;
    }
    else
    {
        return 1. / (2. * epsilon) * (1 + std::cos(Math::pi * x / epsilon));
    }
}

/** Signed minimum
 *
 * Smallest element with the same sign. If both have a different sign, take
 * the signed epsilon.
 */
template <typename T>
T signed_minimum(T a, T b, T epsilon, T sign)
{
    if (sign > 0)
    {
        if ( (a > 0) && (b > 0) )
        {
            return std::min(a, b);
        }
        else if (a > 0)
        {
            return a;
        }
        else if (b > 0)
        {
            return b;
        }
        else
        {
            return epsilon;
        }
    }
    else
    {
        if ( (a < 0) && (b < 0) )
        {
            return std::max(a, b);
        }
        else if (a < 0)
        {
            return a;
        }
        else if (b < 0)
        {
            return b;
        }
        else
        {
            return -1*epsilon;
        }
    }
}

} // namespace itk

#endif // itkSignalFunctions_h
