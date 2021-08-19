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
#ifndef itkFiniteDifference_h
#define itkFiniteDifference_h

#include "itkMath.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkCovariantVector.h"

namespace itk
{

/** Fourth-order accurate, central difference, first derivative
 * The following values are used:
 *  a = f(x - 2 h)
 *  b = f(x - 1 h)
 *  c = f(x + 0 h)
 *  d = f(x + 1 h)
 *  e = f(x + 2 h)
 *  h = spacing
 */
template <typename T>
T first_derivative_central_fourth_order(T a, T b, T c, T d, T e, T h)
{
    T num = 1./12.*a + -2./3.*b + 0.*c + 2./3.*d + -1./12.*e;
    return static_cast<T>(num / h);
}

/** Fourth-order accurate, central difference, second derivative
 * The following values are used:
 *  a = f(x - 2 h)
 *  b = f(x - 1 h)
 *  c = f(x + 0 h)
 *  d = f(x + 1 h)
 *  e = f(x + 2 h)
 *  h = spacing
*/
template <typename T>
T second_derivative_central_fourth_order(T a, T b, T c, T d, T e, T h)
{
    T num = -1./12.*a + 4./3.*b + -5./2.*c + 4./3.*d + -1./12.*e;
    return static_cast<T>(num / (h*h));
}

/** Fourth-order, central difference, mixed second derivative
 * See http://www.holoborodko.com/pavel/2014/11/04/computing-mixed-derivatives-by-finite-differences/
 * The following values are used:
 *  (a, b, c, d)_1 = f( 1, -2), f( 2, -1), f(-2,  1), f(-1,  2)
 *  (a, b, c, d)_2 = f(-1, -2), f(-2, -1), f( 1,  2), f( 2,  1)
 *  (a, b, c, d)_3 = f( 2, -2), f(-2,  2), f(-2, -2), f( 2,  2)
 *  (a, b, c, d)_4 = f(-1, -1), f( 1,  1), f( 1, -1), f(-1,  1)
 *  h_x, h_y = spacing in x and y
 */
template <typename T>
T mixed_second_derivative_central_fourth_order(
    T a1, T b1, T c1, T d1,
    T a2, T b2, T c2, T d2,
    T a3, T b3, T c3, T d3,
    T a4, T b4, T c4, T d4,
    T h_x, T h_y
){
    T num = 0.;
    num += 8.0 * (a1 + b1 + c1 + d1);
    num += -8.0 * (a2 + b2 + c2 + d2);
    num += -1.0 * (a3 + b3 - c3 - d3);
    num += 64.0 * (a4 + b4 - c4 - d4);
    return static_cast<T>(num / (144. * h_x * h_y));
}

/** WENO flux
 * Jiang, Guang-Shan, and Danping Peng.
 * "Weighted ENO schemes for Hamilton--Jacobi equations."
 * SIAM Journal on Scientific computing 21.6 (2000): 2126-2143.
 */
template <typename T>
T weno_fifth(T a, T b, T c, T d, T epsilon)
{

    /* Stencil smoothness */
    T IS0 = 13.0*std::pow(a-b, 2.) + 3.0*std::pow(a-3.0*b, 2.);
    T IS1 = 13.0*std::pow(b-c, 2.) + 3.0*std::pow(b+c, 2.);
    T IS2 = 13.0*std::pow(c-d, 2.) + 3.0*std::pow(3.0*c-d, 2.);

    /* alpha */
    T a0 = 1.0 / std::pow(epsilon + IS0, 2.);
    T a1 = 6.0 / std::pow(epsilon + IS1, 2.);
    T a2 = 3.0 / std::pow(epsilon + IS2, 2.);

    /* weights */
    T w0 = a0 / (a0 + a1 + a2);
    T w2 = a2 / (a0 + a1 + a2);

    return 1.0/3.0*w0*(a-2.0*b+c) + 1.0/6.0*(w2 - 0.5)*(b-2.0*c+d);
}

/** WENO scheme for left-sided first derivative
 * epsilon = 1e-6 is used in paper.
 *
 * Zhang, Yong-Tao, Hong-Kai Zhao, and Jianliang Qian.
 * "High order fast sweeping methods for static Hamilton–Jacobi
 * equations." Journal of Scientific Computing 29.1 (2006): 25-56.
 */
template <typename T>
T weno_negative_third_order(T phi_n_2, T phi_n_1, T phi_0, T phi_p_1, T h, T epsilon)
{
    /* w- */
    T r_num = epsilon + std::pow(phi_0 - 2.*phi_n_1 + phi_n_2, 2);
    T r_den = epsilon + std::pow(phi_p_1 - 2.*phi_0 + phi_n_1, 2);
    T r = r_num / r_den;
    T w = 1. / (1. + 2. * std::pow(r, 2));

    /* Sided difference */
    T first = (phi_p_1 - phi_n_1) / (2.*h);
    T second = (3.*phi_0 - 4.*phi_n_1 + phi_n_2) / (2.*h);
    return (1. - w) * first + w * second;
}

/** WENO scheme for right-sided first derivative
 * epsilon = 1e-6 is used in paper.
 *
 * Zhang, Yong-Tao, Hong-Kai Zhao, and Jianliang Qian.
 * "High order fast sweeping methods for static Hamilton–Jacobi
 * equations." Journal of Scientific Computing 29.1 (2006): 25-56.
 */
template <typename T>
T weno_positive_third_order(T phi_n_1, T phi_0, T phi_p_1, T phi_p_2, T h, T epsilon)
{
    /* w+ */
    T r_num = epsilon + std::pow(phi_p_2 - 2.*phi_p_1 + phi_0, 2);
    T r_den = epsilon + std::pow(phi_p_1 - 2.*phi_0 + phi_n_1, 2);
    T r = r_num / r_den;
    T w = 1. / (1. + 2. * std::pow(r, 2));

    /* Sided difference */
    T first = (phi_p_1 - phi_n_1) / (2.*h);
    T second = (-1.*phi_p_2 + 4.*phi_p_1 + -3.*phi_0) / (2.*h);
    return (1. - w) * first + w * second;
}

/** WENO scheme for left-sided first derivative
 * epsilon = 1e-6 is used in paper.
 *
 * Jiang, Guang-Shan, and Danping Peng.
 * "Weighted ENO schemes for Hamilton--Jacobi equations."
 * SIAM Journal on Scientific computing 21.6 (2000): 2126-2143.
 */
template <typename T>
T weno_negative_fifth_order(T phi_n_3, T phi_n_2, T phi_n_1, T phi_0, T phi_p_1, T phi_p_2, T h, T epsilon)
{
    /* Derivatives */
    T d_phi_n_3 = (phi_n_2-phi_n_3)/h;
    T d_phi_n_2 = (phi_n_1-phi_n_2)/h;
    T d_phi_n_1 = (phi_0-phi_n_1)/h;
    T d_phi_0   = (phi_p_1-phi_0)/h;
    T d_phi_p_1 = (phi_p_2-phi_p_1)/h;

    /* Smoothness */
    T a = d_phi_n_2 - d_phi_n_3;
    T b = d_phi_n_1 - d_phi_n_2;
    T c = d_phi_0 - d_phi_n_1;
    T d = d_phi_p_1 - d_phi_0;

    /* Sided difference */
    T derivative = 1.0/12.0*(-1.0*d_phi_n_2 + 7.0*d_phi_n_1 + 7.0*d_phi_0 - 1.0*d_phi_p_1);
    return derivative - weno_fifth<T>(a, b, c, d, epsilon);
}

/** WENO scheme for right-sided first derivative
 * epsilon = 1e-6 is used in paper.
 *
 * Jiang, Guang-Shan, and Danping Peng.
 * "Weighted ENO schemes for Hamilton--Jacobi equations."
 * SIAM Journal on Scientific computing 21.6 (2000): 2126-2143.
 */
template <typename T>
T weno_positive_fifth_order(T phi_n_2, T phi_n_1, T phi_0, T phi_p_1, T phi_p_2, T phi_p_3, T h, T epsilon)
{
    /* Derivatives */
    T d_phi_n_2 = (phi_n_1-phi_n_2)/h;
    T d_phi_n_1 = (phi_0-phi_n_1)/h;
    T d_phi_0   = (phi_p_1-phi_0)/h;
    T d_phi_p_1 = (phi_p_2-phi_p_1)/h;
    T d_phi_p_2 = (phi_p_3-phi_p_2)/h;

    /* Smoothness */
    T a = d_phi_p_2 - d_phi_p_1;
    T b = d_phi_p_1 - d_phi_0;
    T c = d_phi_0 - d_phi_n_1;
    T d = d_phi_n_1 - d_phi_n_2;

    /* Sided difference */
    T derivative = 1.0/12.0*(-1.0*d_phi_n_2 + 7.0*d_phi_n_1 + 7.0*d_phi_0 - 1.0*d_phi_p_1);
    return derivative + weno_fifth<T>(a, b, c, d, epsilon);
}

/** Compute mean curvature from gradient and Hessian
 * An epsilon is available to make the computation stable in division. One
 * can typically take epsilon to be 1e-9.
 *
 * This function divides out the dimensionality of the surface (NDimension - 1)
 * so the mean curvature agrees with the differential geometry definition,
 * not the physics definition.
 */
template<typename T, unsigned int NDimension = 3>
T mean_curvature(CovariantVector<T,NDimension> first_derivatives, SymmetricSecondRankTensor<T,NDimension> second_derivatives, T epsilon)
{
    unsigned int d = NDimension;
    T mean_curvature = 0.;

    /* Compute numerator */
    for (unsigned int i = 0; i < d; i++)
    {
        T second = 0.;
        for (unsigned int j = 0; j < d; j++)
        {
            if (i == j) {continue;}
            second += second_derivatives(j, j);
        }
        mean_curvature += std::pow(first_derivatives[i], 2.) * second;
    }

    for (unsigned int i = 0; i < d; i++)
    {
        for (unsigned int j = 0; j < d; j++)
        {
            mean_curvature -= 2.*first_derivatives[i]*first_derivatives[j]*second_derivatives(j, i);
        }
    }

    /* Compute denominator */
    T mag_grad = 0;
    for (unsigned int i = 0; i < d; i++)
    {
        mag_grad += std::pow(first_derivatives[i], 2.);
    }
    mag_grad = std::pow(mag_grad, 3./2.);

    /* Divide out */
    mean_curvature /= (mag_grad + epsilon);
    mean_curvature /= std::max(1.0, static_cast<T>(d-1));

    return mean_curvature;
}

/** Compute the Gaussian curvature given the gradient and Hessian
 * An epsilon is available to make the computation stable in division. One
 * can typically take epsilon to be 1e-9.
 */
template<typename T, unsigned int NDimension = 3>
T gaussian_curvature(CovariantVector<T,NDimension> first_derivatives, SymmetricSecondRankTensor<T,NDimension> second_derivatives, T epsilon)
{
    unsigned int d = NDimension;
    T gaussian_curvature = 0.;

    /* Compute numerator */
    for (unsigned int i = 0; i < d; i++)
    {
        unsigned int j = (i+1) % d;
        unsigned int k = (i+2) % d;

        gaussian_curvature += std::pow(first_derivatives[i], 2.) * (
            second_derivatives(j, j)*second_derivatives(k, k) - std::pow(second_derivatives(j, k), 2.)
        );
    }
    for (unsigned int i = 0; i < d; i++)
    {
        unsigned int j = (i+1) % d;
        unsigned int k = (i+2) % d;

        gaussian_curvature += 2.*(
            first_derivatives[i]*first_derivatives[j]*(
                second_derivatives(i, k)*second_derivatives(j, k) -
                second_derivatives(i, j)*second_derivatives(k, k)
            )
        );
    }

    /* Compute denominator */
    T mag_grad = 0;
    for (unsigned int i = 0; i < d; i++)
    {
        mag_grad += std::pow(first_derivatives[i], 2.);
    }
    mag_grad = std::pow(mag_grad, 2.);

    /* Divide out */
    gaussian_curvature /= (mag_grad + epsilon);

    return gaussian_curvature;
}

} // namespace itk

#endif // itkFiniteDifference_h
