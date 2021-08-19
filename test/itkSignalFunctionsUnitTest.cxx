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
#include "itkSignalFunctions.h"

#include "itkGTest.h"
#include "itkTestingMacros.h"
#include "itkMath.h"

using TestType = double;

// heaviside_tanh
TEST(heaviside_tanh, TestZero) {
    ASSERT_NEAR(itk::heaviside_tanh<TestType>(0., 1.), 0.5, 1e-9);
}

TEST(heaviside_tanh, TestLargePositiveValue) {
    ASSERT_NEAR(itk::heaviside_tanh<TestType>(1000000., 1.), 1.0, 1e-9);
}

TEST(heaviside_tanh, TestLargeNegativeValue) {
    ASSERT_NEAR(itk::heaviside_tanh<TestType>(-1000000., 1.), 0.0, 1e-9);
}

TEST(heaviside_tanh, TestSpecific) {
    ASSERT_NEAR(itk::heaviside_tanh<TestType>(0.1, 1.), 0.54983399731247795, 1e-9);
}

TEST(heaviside_tanh, TestSpecificEpsilon) {
    ASSERT_NEAR(itk::heaviside_tanh<TestType>(0.1, 0.1), 0.88079707797788243, 1e-9);
}

// dirac_delta_tanh
TEST(dirac_delta_tanh, TestZero) {
    ASSERT_NEAR(itk::dirac_delta_tanh<TestType>(0., 1.), 0.5, 1e-9);
}

TEST(dirac_delta_tanh, TestLargePositiveValue) {
    ASSERT_NEAR(itk::dirac_delta_tanh<TestType>(1000000., 1.), 0.0, 1e-9);
}

TEST(dirac_delta_tanh, TestLargeNegativeValue) {
    ASSERT_NEAR(itk::dirac_delta_tanh<TestType>(-1000000., 1.), 0.0, 1e-9);
}

TEST(dirac_delta_tanh, TestSpecific) {
    ASSERT_NEAR(itk::dirac_delta_tanh<TestType>(0.1, 1.), 0.49503314542371996, 1e-9);
}

TEST(dirac_delta_tanh, TestSpecificEpsilon) {
    ASSERT_NEAR(itk::dirac_delta_tanh<TestType>(0.1, 0.1), 2.0998717080701303, 1e-9);
}

// heaviside_sin
TEST(heaviside_sin, TestZero) {
    ASSERT_NEAR(itk::heaviside_sin<TestType>(0., 1.), 0.5, 1e-9);
}

TEST(heaviside_sin, TestLargePositiveValue) {
    ASSERT_NEAR(itk::heaviside_sin<TestType>(1000000., 1.), 1.0, 1e-9);
}

TEST(heaviside_sin, TestLargeNegativeValue) {
    ASSERT_NEAR(itk::heaviside_sin<TestType>(-1000000., 1.), 0.0, 1e-9);
}

TEST(heaviside_sin, TestSpecific) {
    ASSERT_NEAR(itk::heaviside_sin<TestType>(0.1, 1.), 0.59918158215417339, 1e-9);
}

TEST(heaviside_sin, TestSpecificEpsilon) {
    ASSERT_NEAR(itk::heaviside_sin<TestType>(0.1, 0.1), 1., 1e-9);
}

// dirac_delta_sin
TEST(dirac_delta_sin, TestZero) {
    ASSERT_NEAR(itk::dirac_delta_sin<TestType>(0., 1.), 1.0, 1e-9);
}

TEST(dirac_delta_sin, TestLargePositiveValue) {
    ASSERT_NEAR(itk::dirac_delta_sin<TestType>(1000000., 1.), 0.0, 1e-9);
}

TEST(dirac_delta_sin, TestLargeNegativeValue) {
    ASSERT_NEAR(itk::dirac_delta_sin<TestType>(-1000000., 1.), 0.0, 1e-9);
}

TEST(dirac_delta_sin, TestSpecific) {
    ASSERT_NEAR(itk::dirac_delta_sin<TestType>(0.1, 1.), 0.97552825814757682, 1e-9);
}

TEST(dirac_delta_sin, TestSpecificEpsilon) {
    ASSERT_NEAR(itk::dirac_delta_sin<TestType>(0.1, 0.1), 0., 1e-9);
}

// signed_minimum
TEST(signed_minimum, TestPosPosPos) {
    ASSERT_NEAR(itk::signed_minimum<TestType>(1., 1., 1e-6, +1), 1.0, 1e-9);
}

TEST(signed_minimum, TestPosNegPos) {
    ASSERT_NEAR(itk::signed_minimum<TestType>(1., -1., 1e-6, +1), 1.0, 1e-9);
}

TEST(signed_minimum, TestNegPosPos) {
    ASSERT_NEAR(itk::signed_minimum<TestType>(-1., 1., 1e-6, +1), 1.0, 1e-9);
}

TEST(signed_minimum, TestNegNegPos) {
    ASSERT_NEAR(itk::signed_minimum<TestType>(-1., -1., 1e-6, +1), 1e-6, 1e-9);
}

TEST(signed_minimum, TestNegNegNeg) {
    ASSERT_NEAR(itk::signed_minimum<TestType>(-1., -1., 1e-6, -1), -1.0, 1e-9);
}

TEST(signed_minimum, TestNegPosNeg) {
    ASSERT_NEAR(itk::signed_minimum<TestType>(-1., +1., 1e-6, -1), -1.0, 1e-9);
}

TEST(signed_minimum, TestPosNegNeg) {
    ASSERT_NEAR(itk::signed_minimum<TestType>(+1., -1., 1e-6, -1), -1.0, 1e-9);
}

TEST(signed_minimum, TestPosPosNeg) {
    ASSERT_NEAR(itk::signed_minimum<TestType>(+1., +1., 1e-6, -1), -1e-6, 1e-9);
}
