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
#include "itkFiniteDifference.h"

#include "itkGTest.h"
#include "itkTestingMacros.h"
#include "itkMath.h"

using TestType = double;

// first_derivative_central_fourth_order
TEST(first_derivative_central_fourth_order, TestAllZeros) {
    ASSERT_NEAR(itk::first_derivative_central_fourth_order<TestType>(0., 0., 0., 0., 0., 1.), 0., 1e-9);
}

TEST(first_derivative_central_fourth_order, TestAllOnes) {
    ASSERT_NEAR(itk::first_derivative_central_fourth_order<TestType>(1., 1., 1., 1., 1., 1.), 0., 1e-9);
}

TEST(first_derivative_central_fourth_order, TestPositive) {
    ASSERT_NEAR(itk::first_derivative_central_fourth_order<TestType>(1., 2., 3., 4., 5., 1.), 1., 1e-9);
}

TEST(first_derivative_central_fourth_order, TestPositiveSpacing) {
    ASSERT_NEAR(itk::first_derivative_central_fourth_order<TestType>(1., 2., 3., 4., 5., 0.5), 2., 1e-9);
}

TEST(first_derivative_central_fourth_order, TestNegative) {
    ASSERT_NEAR(itk::first_derivative_central_fourth_order<TestType>(5., 4., 3., 2., 1, 1.0), -1., 1e-9);
}

TEST(first_derivative_central_fourth_order, TestNegativeSpacing) {
    ASSERT_NEAR(itk::first_derivative_central_fourth_order<TestType>(5., 4., 3., 2., 1, 0.5), -2., 1e-9);
}

TEST(first_derivative_central_fourth_order, TestSpecific) {
    ASSERT_NEAR(itk::first_derivative_central_fourth_order<TestType>(1., 0.99500417, 0.98006658, 0.95533649, 0.92106099, 1.0), -0.019866869166666606, 1e-9);
}

TEST(first_derivative_central_fourth_order, TestSpecificSpacing) {
    ASSERT_NEAR(itk::first_derivative_central_fourth_order<TestType>(1., 0.99500417, 0.98006658, 0.95533649, 0.92106099, 2.5), -0.007946747666666643, 1e-9);
}

// second_derivative_central_fourth_order
TEST(second_derivative_central_fourth_order, TestAllZeros) {
    ASSERT_NEAR(itk::second_derivative_central_fourth_order<TestType>(0., 0., 0., 0., 0., 1.), 0., 1e-9);
}

TEST(second_derivative_central_fourth_order, TestAllOnes) {
    ASSERT_NEAR(itk::second_derivative_central_fourth_order<TestType>(1., 1., 1., 1., 1., 1.), 0., 1e-9);
}

TEST(second_derivative_central_fourth_order, TestPositive) {
    ASSERT_NEAR(itk::second_derivative_central_fourth_order<TestType>(1., 4., 9., 16., 25., 1.), 2., 1e-9);
}

TEST(second_derivative_central_fourth_order, TestPositiveSpacing) {
    ASSERT_NEAR(itk::second_derivative_central_fourth_order<TestType>(1., 4., 9., 16., 25., 2.0), 0.5, 1e-9);
}

TEST(second_derivative_central_fourth_order, TestNegative) {
    ASSERT_NEAR(itk::second_derivative_central_fourth_order<TestType>(1.0, 2.0, 4.0, 2.0, 1.0, 1.0), -4.833333333333333, 1e-9);
}

TEST(second_derivative_central_fourth_order, TestNegativeSpacing) {
    ASSERT_NEAR(itk::second_derivative_central_fourth_order<TestType>(1.0, 2.0, 4.0, 2.0, 1.0, 2.0), -1.2083333333333333, 1e-9);
}

TEST(second_derivative_central_fourth_order, TestSpecific) {
    ASSERT_NEAR(itk::second_derivative_central_fourth_order<TestType>(1., 0.99500417, 0.98006658, 0.95533649, 0.92106099, 1.0), -0.009800652500000201, 1e-9);
}

TEST(second_derivative_central_fourth_order, TestSpecificSpacing) {
    ASSERT_NEAR(itk::second_derivative_central_fourth_order<TestType>(1., 0.99500417, 0.98006658, 0.95533649, 0.92106099, 2.5), -0.001568104400000032, 1e-9);
}

// mixed_second_derivative_central_fourth_order
TEST(mixed_second_derivative_central_fourth_order, TestAllZeros) {
    TestType result = itk::mixed_second_derivative_central_fourth_order<TestType>(
      0., 0., 0., 0.,
      0., 0., 0., 0.,
      0., 0., 0., 0.,
      0., 0., 0., 0.,
      1., 1.
    );
    ASSERT_NEAR(result, 0., 1e-9);
}

TEST(mixed_second_derivative_central_fourth_order, TestAllOnes) {
    TestType result = itk::mixed_second_derivative_central_fourth_order<TestType>(
      1., 1., 1., 1.,
      1., 1., 1., 1.,
      1., 1., 1., 1.,
      1., 1., 1., 1.,
      1., 1.
    );
    ASSERT_NEAR(result, 0., 1e-9);
}

TEST(mixed_second_derivative_central_fourth_order, TestPositive) {
    TestType result = itk::mixed_second_derivative_central_fourth_order<TestType>(
      1., 2., 3., 4.,
      2., 2., 2., 2.,
      4., 4., 4., 4.,
      1., 1., 1., 1.,
      1., 1.
    );
    ASSERT_NEAR(result, 0.1111111111111111, 1e-9);
}

TEST(mixed_second_derivative_central_fourth_order, TestPositiveSpacingX) {
    TestType result = itk::mixed_second_derivative_central_fourth_order<TestType>(
      1., 2., 3., 4.,
      2., 2., 2., 2.,
      4., 4., 4., 4.,
      1., 1., 1., 1.,
      2., 1.
    );
    ASSERT_NEAR(result, 0.055555555555555552, 1e-9);
}

TEST(mixed_second_derivative_central_fourth_order, TestPositiveSpacingY) {
    TestType result = itk::mixed_second_derivative_central_fourth_order<TestType>(
      1., 2., 3., 4.,
      2., 2., 2., 2.,
      4., 4., 4., 4.,
      1., 1., 1., 1.,
      1., 2.
    );
    ASSERT_NEAR(result, 0.055555555555555552, 1e-9);
}

TEST(mixed_second_derivative_central_fourth_order, TestPositiveSpacing) {
    TestType result = itk::mixed_second_derivative_central_fourth_order<TestType>(
      1., 2., 3., 4.,
      2., 2., 2., 2.,
      4., 4., 4., 4.,
      1., 1., 1., 1.,
      2., 2.
    );
    ASSERT_NEAR(result, 0.027777777777777776, 1e-9);
}

TEST(mixed_second_derivative_central_fourth_order, TestNegative) {
    TestType result = itk::mixed_second_derivative_central_fourth_order<TestType>(
      -1., -2., -3., -4.,
      -2., -2., -2., -2.,
      -4., -4., -4., -4.,
      -1., -1., -1., -1.,
      1., 1.
    );
    ASSERT_NEAR(result, -0.1111111111111111, 1e-9);
}

TEST(mixed_second_derivative_central_fourth_order, TestNegativeSpacingX) {
    TestType result = itk::mixed_second_derivative_central_fourth_order<TestType>(
      -1., -2., -3., -4.,
      -2., -2., -2., -2.,
      -4., -4., -4., -4.,
      -1., -1., -1., -1.,
      2., 1.
    );
    ASSERT_NEAR(result, -0.055555555555555552, 1e-9);
}

TEST(mixed_second_derivative_central_fourth_order, TestNegativeSpacingY) {
    TestType result = itk::mixed_second_derivative_central_fourth_order<TestType>(
      -1., -2., -3., -4.,
      -2., -2., -2., -2.,
      -4., -4., -4., -4.,
      -1., -1., -1., -1.,
      1., 2.
    );
    ASSERT_NEAR(result, -0.055555555555555552, 1e-9);
}

TEST(mixed_second_derivative_central_fourth_order, TestNegativeSpacing) {
    TestType result = itk::mixed_second_derivative_central_fourth_order<TestType>(
      -1., -2., -3., -4.,
      -2., -2., -2., -2.,
      -4., -4., -4., -4.,
      -1., -1., -1., -1.,
      2., 2.
    );
    ASSERT_NEAR(result, -0.027777777777777776, 1e-9);
}

TEST(mixed_second_derivative_central_fourth_order, TestSpecific) {
    TestType result = itk::mixed_second_derivative_central_fourth_order<TestType>(
      0.55646719, 0.50018713, 0.13374574, 0.04221229,
      0.0052251 , 0.44507512, 0.46862878, 0.13483009,
      0.10676733, 0.97073443, 0.83308078, 0.53278424,
      0.15779718, 0.4785072 , 0.34206369, 0.40353729,
      1., 1.
    );
    ASSERT_NEAR(result, -0.036637451805555547, 1e-9);
}

TEST(mixed_second_derivative_central_fourth_order, TestSpecificSpacing) {
    TestType result = itk::mixed_second_derivative_central_fourth_order<TestType>(
      0.55646719, 0.50018713, 0.13374574, 0.04221229,
      0.0052251 , 0.44507512, 0.46862878, 0.13483009,
      0.10676733, 0.97073443, 0.83308078, 0.53278424,
      0.15779718, 0.4785072 , 0.34206369, 0.40353729,
      3., 4.
    );
    ASSERT_NEAR(result, -0.0030531209837962956, 1e-9);
}

TEST(mixed_second_derivative_central_fourth_order, TestSpecificSpacingX) {
    TestType result = itk::mixed_second_derivative_central_fourth_order<TestType>(
      0.55646719, 0.50018713, 0.13374574, 0.04221229,
      0.0052251 , 0.44507512, 0.46862878, 0.13483009,
      0.10676733, 0.97073443, 0.83308078, 0.53278424,
      0.15779718, 0.4785072 , 0.34206369, 0.40353729,
      3., 1.
    );
    ASSERT_NEAR(result, -0.012212483935185182, 1e-9);
}

TEST(mixed_second_derivative_central_fourth_order, TestSpecificSpacingY) {
    TestType result = itk::mixed_second_derivative_central_fourth_order<TestType>(
      0.55646719, 0.50018713, 0.13374574, 0.04221229,
      0.0052251 , 0.44507512, 0.46862878, 0.13483009,
      0.10676733, 0.97073443, 0.83308078, 0.53278424,
      0.15779718, 0.4785072 , 0.34206369, 0.40353729,
      1., 4.
    );
    ASSERT_NEAR(result, -0.0091593629513888868, 1e-9);
}

// weno_negative_fifth_order
TEST(weno_negative_fifth_order, TestAllZeros) {
    ASSERT_NEAR(itk::weno_negative_fifth_order<TestType>(0., 0., 0., 0., 0., 0., 1., 1e-6), 0., 1e-9);
}

TEST(weno_negative_fifth_order, TestAllOnes) {
    ASSERT_NEAR(itk::weno_negative_fifth_order<TestType>(1., 1., 1., 1., 1., 0., 1., 1e-6), 0., 1e-9);
}

TEST(weno_negative_fifth_order, TestPositive) {
    ASSERT_NEAR(itk::weno_negative_fifth_order<TestType>(1., 2., 3., 4., 5., 6., 1., 1e-6), 1., 1e-9);
}

TEST(weno_negative_fifth_order, TestPositiveSpacing) {
    ASSERT_NEAR(itk::weno_negative_fifth_order<TestType>(1., 2., 3., 4., 5., 6., 0.5, 1e-6), 2., 1e-9);
}

TEST(weno_negative_fifth_order, TestNegative) {
    ASSERT_NEAR(itk::weno_negative_fifth_order<TestType>(5., 4., 3., 2., 1, 0., 1.0, 1e-6), -1., 1e-9);
}

TEST(weno_negative_fifth_order, TestNegativeSpacing) {
    ASSERT_NEAR(itk::weno_negative_fifth_order<TestType>(5., 4., 3., 2., 1, 0., 0.5, 1e-6), -2., 1e-9);
}

TEST(weno_negative_fifth_order, TestSpecific) {
    ASSERT_NEAR(itk::weno_negative_fifth_order<TestType>(1., 0.99500417, 0.98006658, 0.95533649, 0.92106099, 0.88123432123, 1.0, 1e-6), -0.029636883653833566, 1e-9);
}

TEST(weno_negative_fifth_order, TestSpecificSpacing) {
    ASSERT_NEAR(itk::weno_negative_fifth_order<TestType>(1., 0.99500417, 0.98006658, 0.95533649, 0.92106099, 0.88123432123, 2.5, 1e-6), -0.01185486000703519, 1e-9);
}

// weno_positive_fifth_order
TEST(weno_positive_fifth_order, TestAllZeros) {
    ASSERT_NEAR(itk::weno_positive_fifth_order<TestType>(0., 0., 0., 0., 0., 0., 1., 1e-6), 0., 1e-9);
}

TEST(weno_positive_fifth_order, TestAllOnes) {
    ASSERT_NEAR(itk::weno_positive_fifth_order<TestType>(1., 1., 1., 1., 1., 0., 1., 1e-6), 0., 1e-9);
}

TEST(weno_positive_fifth_order, TestPositive) {
    ASSERT_NEAR(itk::weno_positive_fifth_order<TestType>(1., 2., 3., 4., 5., 6., 1., 1e-6), 1., 1e-9);
}

TEST(weno_positive_fifth_order, TestPositiveSpacing) {
    ASSERT_NEAR(itk::weno_positive_fifth_order<TestType>(1., 2., 3., 4., 5., 6., 0.5, 1e-6), 2., 1e-9);
}

TEST(weno_positive_fifth_order, TestNegative) {
    ASSERT_NEAR(itk::weno_positive_fifth_order<TestType>(5., 4., 3., 2., 1, 0., 1.0, 1e-6), -1., 1e-9);
}

TEST(weno_positive_fifth_order, TestNegativeSpacing) {
    ASSERT_NEAR(itk::weno_positive_fifth_order<TestType>(5., 4., 3., 2., 1, 0., 0.5, 1e-6), -2., 1e-9);
}

TEST(weno_positive_fifth_order, TestSpecific) {
    ASSERT_NEAR(itk::weno_positive_fifth_order<TestType>(1., 0.99500417, 0.98006658, 0.95533649, 0.92106099, 0.88123432123, 1.0, 1e-6), -0.019818117680131954, 1e-9);
}

TEST(weno_positive_fifth_order, TestSpecificSpacing) {
    ASSERT_NEAR(itk::weno_positive_fifth_order<TestType>(1., 0.99500417, 0.98006658, 0.95533649, 0.92106099, 0.88123432123, 2.5, 1e-6), -0.0079271778408092543, 1e-9);
}

// mean_curvature
TEST(mean_curvature, TestIdentity) {
    const unsigned int NDimension = 3;
    using VectorType = itk::CovariantVector<TestType, NDimension>;
    using TensorType = itk::SymmetricSecondRankTensor<TestType, NDimension>;

    VectorType first_der = VectorType();
    first_der[0] = 1.;
    first_der[1] = 1.;
    first_der[2] = 1.;

    TensorType hessian = TensorType();
    hessian(0, 0) = 1.;
    hessian(0, 1) = 0.;
    hessian(0, 2) = 0.;
    hessian(1, 1) = 1.;
    hessian(1, 2) = 0.;
    hessian(2, 2) = 1.;

    ASSERT_NEAR(itk::mean_curvature<TestType>(first_der, hessian, 1e-9), 0., 1e-9);
}

TEST(mean_curvature, TestZeros) {
    const unsigned int NDimension = 3;
    using VectorType = itk::CovariantVector<TestType, NDimension>;
    using TensorType = itk::SymmetricSecondRankTensor<TestType, NDimension>;

    VectorType first_der = VectorType();
    first_der[0] = 0.;
    first_der[1] = 0.;
    first_der[2] = 0.;

    TensorType hessian = TensorType();
    hessian(0, 0) = 0.;
    hessian(0, 1) = 0.;
    hessian(0, 2) = 0.;
    hessian(1, 1) = 0.;
    hessian(1, 2) = 0.;
    hessian(2, 2) = 0.;

    ASSERT_NEAR(itk::mean_curvature<TestType>(first_der, hessian, 1e-9), 0., 1e-9);
}

TEST(mean_curvature, TestSpecific) {
    const unsigned int NDimension = 3;
    using VectorType = itk::CovariantVector<TestType, NDimension>;
    using TensorType = itk::SymmetricSecondRankTensor<TestType, NDimension>;

    VectorType first_der = VectorType();
    first_der[0] = 1.;
    first_der[1] = 2.;
    first_der[2] = 3.;

    TensorType hessian = TensorType();
    hessian(0, 0) = 0.;
    hessian(0, 1) = 0.;
    hessian(0, 2) = 0.;
    hessian(1, 1) = 0.1542515;
    hessian(1, 2) = 0.9880316;
    hessian(2, 2) = 0.1542515;

    ASSERT_NEAR(itk::mean_curvature<TestType>(first_der, hessian, 1e-9), -0.24253504218032254, 1e-9);
}

// gaussian_curvature
TEST(gaussian_curvature, TestIdentity) {
    const unsigned int NDimension = 3;
    using VectorType = itk::CovariantVector<TestType, NDimension>;
    using TensorType = itk::SymmetricSecondRankTensor<TestType, NDimension>;

    VectorType first_der = VectorType();
    first_der[0] = 1.;
    first_der[1] = 1.;
    first_der[2] = 1.;

    TensorType hessian = TensorType();
    hessian(0, 0) = 1.;
    hessian(0, 1) = 0.;
    hessian(0, 2) = 0.;
    hessian(1, 1) = 1.;
    hessian(1, 2) = 0.;
    hessian(2, 2) = 1.;

    ASSERT_NEAR(itk::gaussian_curvature<TestType>(first_der, hessian, 1e-9), 0.33333333329629627, 1e-9);
}

TEST(gaussian_curvature, TestZeros) {
    const unsigned int NDimension = 3;
    using VectorType = itk::CovariantVector<TestType, NDimension>;
    using TensorType = itk::SymmetricSecondRankTensor<TestType, NDimension>;

    VectorType first_der = VectorType();
    first_der[0] = 0.;
    first_der[1] = 0.;
    first_der[2] = 0.;

    TensorType hessian = TensorType();
    hessian(0, 0) = 0.;
    hessian(0, 1) = 0.;
    hessian(0, 2) = 0.;
    hessian(1, 1) = 0.;
    hessian(1, 2) = 0.;
    hessian(2, 2) = 0.;

    ASSERT_NEAR(itk::gaussian_curvature<TestType>(first_der, hessian, 1e-9), 0., 1e-9);
}

TEST(gaussian_curvature, TestSpecific) {
    const unsigned int NDimension = 3;
    using VectorType = itk::CovariantVector<TestType, NDimension>;
    using TensorType = itk::SymmetricSecondRankTensor<TestType, NDimension>;

    VectorType first_der = VectorType();
    first_der[0] = 1.;
    first_der[1] = 2.;
    first_der[2] = 3.;

    TensorType hessian = TensorType();
    hessian(0, 0) = 0.;
    hessian(0, 1) = 0.;
    hessian(0, 2) = 0.;
    hessian(1, 1) = 0.1542515;
    hessian(1, 2) = 0.9880316;
    hessian(2, 2) = 0.1542515;

    ASSERT_NEAR(itk::gaussian_curvature<TestType>(first_der, hessian, 1e-9), -0.0048592495782727087, 1e-9);
}
