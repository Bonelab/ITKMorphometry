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
#include "itkHeavisideImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkGTest.h"
#include "itkTestingMacros.h"
#include "itkMath.h"

namespace
{
template <typename T>
class itkHeavisideImageFilterUnitTest
  : public ::testing::Test
{
public:
  /* Useful typedefs */
  static const unsigned int DIMENSION = 3;
  using PixelType         = T;
  using InputImageType    = itk::Image< PixelType, DIMENSION >;
  using FilterType        = typename itk::HeavisideImageFilter< InputImageType >;
  using FilterPointerType = typename FilterType::Pointer;

  itkHeavisideImageFilterUnitTest() {
    /* Instantiate filter */
    m_Filter = FilterType::New();
  }
  ~itkHeavisideImageFilterUnitTest() override {}

protected:
  void SetUp() override {}
  void TearDown() override {}

  FilterPointerType                   m_Filter;
};
}

// Define the templates we would like to test
using TestingLabelTypes = ::testing::Types<double, float>;

TYPED_TEST_SUITE(itkHeavisideImageFilterUnitTest, TestingLabelTypes);

TYPED_TEST(itkHeavisideImageFilterUnitTest, InitialApproximation) {
    using FilterType  = typename TestFixture::FilterType;
    ASSERT_EQ(this->m_Filter->GetApproximation(), itk::SmoothApproximationType::Tanh);
}

TYPED_TEST(itkHeavisideImageFilterUnitTest, InitialEpsilon) {
    ASSERT_EQ(this->m_Filter->GetEpsilon(), 1.);
}

TYPED_TEST(itkHeavisideImageFilterUnitTest, SetEpsilon) {
    this->m_Filter->SetEpsilon(0.01);
    ASSERT_EQ(this->m_Filter->GetEpsilon(), 0.01);
}

TYPED_TEST(itkHeavisideImageFilterUnitTest, PrintTanh) {
    using FilterType  = typename TestFixture::FilterType;
    ASSERT_EQ(this->m_Filter->GetApproximationAsString(itk::SmoothApproximationType::Tanh), "tanh");
}

TYPED_TEST(itkHeavisideImageFilterUnitTest, PrintSin) {
    using FilterType  = typename TestFixture::FilterType;
    ASSERT_EQ(this->m_Filter->GetApproximationAsString(itk::SmoothApproximationType::Sin), "sin");
}

TYPED_TEST(itkHeavisideImageFilterUnitTest, SetApproximationToTanh) {
    using FilterType  = typename TestFixture::FilterType;
    this->m_Filter->SetApproximationToTanh();
    ASSERT_EQ(this->m_Filter->GetApproximation(), itk::SmoothApproximationType::Tanh);
}

TYPED_TEST(itkHeavisideImageFilterUnitTest, SetApproximationToSin) {
    using FilterType  = typename TestFixture::FilterType;
    this->m_Filter->SetApproximationToSin();
    ASSERT_EQ(this->m_Filter->GetApproximation(), itk::SmoothApproximationType::Sin);
}

TYPED_TEST(itkHeavisideImageFilterUnitTest, ComputeImageTanh) {
    /* Create Image */
    using PixelType  = typename TestFixture::PixelType;
    using InputImageType = typename TestFixture::InputImageType;
    typename InputImageType::Pointer image = InputImageType::New();

    typename InputImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;

    typename InputImageType::SizeType size;
    size[0] = 1;
    size[1] = 1;
    size[2] = 5;

    typename InputImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);

    image->SetRegions(region);
    image->Allocate();
    image->FillBuffer(0.);

    itk::ImageRegionIterator< InputImageType > it(image, image->GetLargestPossibleRegion());

    PixelType f = -1.0;
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        it.Set(f);
        f += 0.5;
    }

    /* Process */
    this->m_Filter->SetEpsilon(1.0);
    this->m_Filter->SetApproximationToTanh();
    this->m_Filter->SetInput(image);
    EXPECT_NO_THROW(this->m_Filter->Update());
    EXPECT_TRUE(this->m_Filter->GetOutput()->GetBufferedRegion() == region);

    /* Check output */
    PixelType expected [5] = { 0.11920291930437088, 0.26894143223762512, 0.5, 0.73105859756469727, 0.88079708814620972 };
    itk::ImageRegionIterator< InputImageType > out(this->m_Filter->GetOutput(), this->m_Filter->GetOutput()->GetLargestPossibleRegion());
    unsigned int i = 0;
    for (out.GoToBegin(); !out.IsAtEnd(); ++out)
    {
        ASSERT_NEAR(out.Get(), expected[i], 1e-6);
        ++i;
    }
    ASSERT_EQ(i, 5);
}

TYPED_TEST(itkHeavisideImageFilterUnitTest, ComputeImageTanhEpsilon) {
    /* Create Image */
    using PixelType  = typename TestFixture::PixelType;
    using InputImageType = typename TestFixture::InputImageType;
    typename InputImageType::Pointer image = InputImageType::New();

    typename InputImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;

    typename InputImageType::SizeType size;
    size[0] = 1;
    size[1] = 1;
    size[2] = 5;

    typename InputImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);

    image->SetRegions(region);
    image->Allocate();
    image->FillBuffer(0.);

    itk::ImageRegionIterator< InputImageType > it(image, image->GetLargestPossibleRegion());

    PixelType f = -1.0;
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        it.Set(f);
        f += 0.5;
    }

    /* Process */
    this->m_Filter->SetEpsilon(2.0);
    this->m_Filter->SetApproximationToTanh();
    this->m_Filter->SetInput(image);
    EXPECT_NO_THROW(this->m_Filter->Update());
    EXPECT_TRUE(this->m_Filter->GetOutput()->GetBufferedRegion() == region);

    /* Check output */
    PixelType expected [5] = { 0.26894143223762512, 0.37754067778587341, 0.5, 0.62245935201644897, 0.73105859756469727 };
    itk::ImageRegionIterator< InputImageType > out(this->m_Filter->GetOutput(), this->m_Filter->GetOutput()->GetLargestPossibleRegion());
    unsigned int i = 0;
    for (out.GoToBegin(); !out.IsAtEnd(); ++out)
    {
        ASSERT_NEAR(out.Get(), expected[i], 1e-6);
        ++i;
    }
    ASSERT_EQ(i, 5);
}

TYPED_TEST(itkHeavisideImageFilterUnitTest, ComputeImageSin) {
    /* Create Image */
    using PixelType  = typename TestFixture::PixelType;
    using InputImageType = typename TestFixture::InputImageType;
    typename InputImageType::Pointer image = InputImageType::New();

    typename InputImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;

    typename InputImageType::SizeType size;
    size[0] = 1;
    size[1] = 1;
    size[2] = 5;

    typename InputImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);

    image->SetRegions(region);
    image->Allocate();
    image->FillBuffer(0.);

    itk::ImageRegionIterator< InputImageType > it(image, image->GetLargestPossibleRegion());

    PixelType f = -1.0;
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        it.Set(f);
        f += 0.5;
    }

    /* Process */
    this->m_Filter->SetEpsilon(1.0);
    this->m_Filter->SetApproximationToSin();
    this->m_Filter->SetInput(image);
    EXPECT_NO_THROW(this->m_Filter->Update());
    EXPECT_TRUE(this->m_Filter->GetOutput()->GetBufferedRegion() == region);

    /* Check output */
    PixelType expected [5] = { 0., 0.090845056908104654, 0.5, 0.90915495157241821, 1. };
    itk::ImageRegionIterator< InputImageType > out(this->m_Filter->GetOutput(), this->m_Filter->GetOutput()->GetLargestPossibleRegion());
    unsigned int i = 0;
    for (out.GoToBegin(); !out.IsAtEnd(); ++out)
    {
        ASSERT_NEAR(out.Get(), expected[i], 1e-6);
        ++i;
    }
    ASSERT_EQ(i, 5);
}

TYPED_TEST(itkHeavisideImageFilterUnitTest, ComputeImageSinEpsilon) {
    /* Create Image */
    using PixelType  = typename TestFixture::PixelType;
    using InputImageType = typename TestFixture::InputImageType;
    typename InputImageType::Pointer image = InputImageType::New();

    typename InputImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;

    typename InputImageType::SizeType size;
    size[0] = 1;
    size[1] = 1;
    size[2] = 5;

    typename InputImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);

    image->SetRegions(region);
    image->Allocate();
    image->FillBuffer(0.);

    itk::ImageRegionIterator< InputImageType > it(image, image->GetLargestPossibleRegion());

    PixelType f = -1.0;
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        it.Set(f);
        f += 0.5;
    }

    /* Process */
    this->m_Filter->SetEpsilon(2.0);
    this->m_Filter->SetApproximationToSin();
    this->m_Filter->SetInput(image);
    EXPECT_NO_THROW(this->m_Filter->Update());
    EXPECT_TRUE(this->m_Filter->GetOutput()->GetBufferedRegion() == region);

    /* Check output */
    PixelType expected [5] = { 0.090845055878162384, 0.26246047019958496, 0.5, 0.73753952980041504, 0.90915495157241821 };
    itk::ImageRegionIterator< InputImageType > out(this->m_Filter->GetOutput(), this->m_Filter->GetOutput()->GetLargestPossibleRegion());
    unsigned int i = 0;
    for (out.GoToBegin(); !out.IsAtEnd(); ++out)
    {
        ASSERT_NEAR(out.Get(), expected[i], 1e-6);
        ++i;
    }
    ASSERT_EQ(i, 5);
}
