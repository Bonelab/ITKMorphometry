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
#include "itkDiracDeltaImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkGTest.h"
#include "itkTestingMacros.h"
#include "itkMath.h"

namespace
{
template <typename T>
class itkDiracDeltaImageFilterUnitTest
  : public ::testing::Test
{
public:
  /* Useful typedefs */
  static const unsigned int DIMENSION = 3;
  using PixelType         = T;
  using InputImageType    = itk::Image< PixelType, DIMENSION >;
  using FilterType        = typename itk::DiracDeltaImageFilter< InputImageType >;
  using FilterPointerType = typename FilterType::Pointer;

  itkDiracDeltaImageFilterUnitTest() {
    /* Instantiate filter */
    m_Filter = FilterType::New();
  }
  ~itkDiracDeltaImageFilterUnitTest() override {}

protected:
  void SetUp() override {}
  void TearDown() override {}

  FilterPointerType                   m_Filter;
};
}

// Define the templates we would like to test
using TestingLabelTypes = ::testing::Types<double, float>;

TYPED_TEST_SUITE(itkDiracDeltaImageFilterUnitTest, TestingLabelTypes);

TYPED_TEST(itkDiracDeltaImageFilterUnitTest, InitialApproximation) {
    using FilterType  = typename TestFixture::FilterType;
    ASSERT_EQ(this->m_Filter->GetApproximation(), FilterType::DiracDeltaApproximationType::Tanh);
}

TYPED_TEST(itkDiracDeltaImageFilterUnitTest, InitialEpsilon) {
    ASSERT_EQ(this->m_Filter->GetEpsilon(), 1.);
}

TYPED_TEST(itkDiracDeltaImageFilterUnitTest, SetEpsilon) {
    this->m_Filter->SetEpsilon(0.01);
    ASSERT_EQ(this->m_Filter->GetEpsilon(), 0.01);
}

TYPED_TEST(itkDiracDeltaImageFilterUnitTest, PrintTanh) {
    using FilterType  = typename TestFixture::FilterType;
    ASSERT_EQ(this->m_Filter->GetApproximationAsString(FilterType::DiracDeltaApproximationType::Tanh), "tanh");
}

TYPED_TEST(itkDiracDeltaImageFilterUnitTest, PrintSin) {
    using FilterType  = typename TestFixture::FilterType;
    ASSERT_EQ(this->m_Filter->GetApproximationAsString(FilterType::DiracDeltaApproximationType::Sin), "sin");
}

TYPED_TEST(itkDiracDeltaImageFilterUnitTest, SetApproximationToTanh) {
    using FilterType  = typename TestFixture::FilterType;
    this->m_Filter->SetApproximationToTanh();
    ASSERT_EQ(this->m_Filter->GetApproximation(), FilterType::DiracDeltaApproximationType::Tanh);
}

TYPED_TEST(itkDiracDeltaImageFilterUnitTest, SetApproximationToSin) {
    using FilterType  = typename TestFixture::FilterType;
    this->m_Filter->SetApproximationToSin();
    ASSERT_EQ(this->m_Filter->GetApproximation(), FilterType::DiracDeltaApproximationType::Sin);
}

TYPED_TEST(itkDiracDeltaImageFilterUnitTest, ComputeImageTanh) {
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
    PixelType expected [5] = { 0.20998716354370117, 0.39322385191917419, 0.5, 0.39322385191917419, 0.20998716354370117 };
    itk::ImageRegionIterator< InputImageType > out(this->m_Filter->GetOutput(), this->m_Filter->GetOutput()->GetLargestPossibleRegion());
    unsigned int i = 0;
    for (out.GoToBegin(); !out.IsAtEnd(); ++out)
    {
        ASSERT_NEAR(out.Get(), expected[i], 1e-6);
        ++i;
    }
    ASSERT_EQ(i, 5);
}

TYPED_TEST(itkDiracDeltaImageFilterUnitTest, ComputeImageTanhEpsilon) {
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
    PixelType expected [5] = { 0.1966119259595871, 0.23500370979309082, 0.25, 0.23500370979309082, 0.1966119259595871 };
    itk::ImageRegionIterator< InputImageType > out(this->m_Filter->GetOutput(), this->m_Filter->GetOutput()->GetLargestPossibleRegion());
    unsigned int i = 0;
    for (out.GoToBegin(); !out.IsAtEnd(); ++out)
    {
        ASSERT_NEAR(out.Get(), expected[i], 1e-6);
        ++i;
    }
    ASSERT_EQ(i, 5);
}

TYPED_TEST(itkDiracDeltaImageFilterUnitTest, ComputeImageSin) {
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
    PixelType expected [5] = { 0., 0.5, 1.0, 0.5, 0.};
    itk::ImageRegionIterator< InputImageType > out(this->m_Filter->GetOutput(), this->m_Filter->GetOutput()->GetLargestPossibleRegion());
    unsigned int i = 0;
    for (out.GoToBegin(); !out.IsAtEnd(); ++out)
    {
        ASSERT_NEAR(out.Get(), expected[i], 1e-6);
        ++i;
    }
    ASSERT_EQ(i, 5);
}

TYPED_TEST(itkDiracDeltaImageFilterUnitTest, ComputeImageSinEpsilon) {
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
    PixelType expected [5] = { 0.25, 0.4267767071723938, 0.5, 0.4267767071723938, 0.25};
    itk::ImageRegionIterator< InputImageType > out(this->m_Filter->GetOutput(), this->m_Filter->GetOutput()->GetLargestPossibleRegion());
    unsigned int i = 0;
    for (out.GoToBegin(); !out.IsAtEnd(); ++out)
    {
        ASSERT_NEAR(out.Get(), expected[i], 1e-6);
        ++i;
    }
    ASSERT_EQ(i, 5);
}
