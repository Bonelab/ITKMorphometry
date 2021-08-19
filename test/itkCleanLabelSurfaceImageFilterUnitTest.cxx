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
#include "itkCleanLabelSurfaceImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkGTest.h"
#include "itkTestingMacros.h"
#include "itkMath.h"

namespace
{
class itkCleanLabelSurfaceImageFilterUnitTest
  : public ::testing::Test
{
public:
  /* Useful typedefs */
  static const unsigned int DIMENSION = 3;
  using PixelType         = unsigned int;
  using InputImageType    = itk::Image< PixelType, DIMENSION >;
  using FilterType        = typename itk::CleanLabelSurfaceImageFilter< InputImageType >;
  using FilterPointerType = typename FilterType::Pointer;

  itkCleanLabelSurfaceImageFilterUnitTest() {
    /* Instantiate filter */
    m_Filter = FilterType::New();
  }
  ~itkCleanLabelSurfaceImageFilterUnitTest() override {}

protected:
  void SetUp() override {}
  void TearDown() override {}

  FilterPointerType                   m_Filter;
};
}

TEST(itkCleanLabelSurfaceImageFilterUnitTest, InitialNumberOfLabels) {
  using PixelType  = unsigned int;
  using InputImageType = itk::Image< PixelType, 3 >;
  using FilterType        = typename itk::CleanLabelSurfaceImageFilter< InputImageType >;
  using FilterPointerType = typename FilterType::Pointer;
  FilterPointerType Filter = FilterType::New();
  ASSERT_EQ(Filter->GetNumberOfLabels(), 0);
}

TEST(itkCleanLabelSurfaceImageFilterUnitTest, InitialBackgroundLabel) {
  using PixelType  = unsigned int;
  using InputImageType = itk::Image< PixelType, 3 >;
  using FilterType        = typename itk::CleanLabelSurfaceImageFilter< InputImageType >;
  using FilterPointerType = typename FilterType::Pointer;
  FilterPointerType Filter = FilterType::New();
  ASSERT_EQ(Filter->GetBackgroundLabel(), 0);
}

TEST(itkCleanLabelSurfaceImageFilterUnitTest, InitialMaxNumberOfLabels) {
  using PixelType  = unsigned int;
  using InputImageType = itk::Image< PixelType, 3 >;
  using FilterType        = typename itk::CleanLabelSurfaceImageFilter< InputImageType >;
  using FilterPointerType = typename FilterType::Pointer;
  FilterPointerType Filter = FilterType::New();
  ASSERT_EQ(Filter->GetMaxNumberOfLabels(), 256);
}

TEST(itkCleanLabelSurfaceImageFilterUnitTest, AllOnesReturnIdentity) {
    /* Create Image */
    using PixelType  = unsigned int;
    using InputImageType = itk::Image< PixelType, 3 >;

    typename InputImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;

    typename InputImageType::SizeType size;
    size[0] = 5;
    size[1] = 5;
    size[2] = 5;

    typename InputImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);

    typename InputImageType::Pointer image = InputImageType::New();
    image->SetRegions(region);
    image->Allocate();
    image->FillBuffer(1);

    /* Filter */
    using FilterType        = typename itk::CleanLabelSurfaceImageFilter< InputImageType >;
    using FilterPointerType = typename FilterType::Pointer;
    FilterPointerType Filter = FilterType::New();
    Filter->SetInput(image);
    Filter->Update();
    EXPECT_NO_THROW(Filter->Update());

    /* Check */
    itk::ImageRegionIterator< InputImageType > input(Filter->GetOutput(), Filter->GetOutput()->GetLargestPossibleRegion());

    for(input.GoToBegin(); !input.IsAtEnd(); ++input)
    {
        ASSERT_EQ(input.Get(), 1);
    }
}

TEST(itkCleanLabelSurfaceImageFilterUnitTest, AllZerosReturnIdentity) {
    /* Create Image */
    using PixelType  = unsigned int;
    using InputImageType = itk::Image< PixelType, 3 >;

    typename InputImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;

    typename InputImageType::SizeType size;
    size[0] = 5;
    size[1] = 5;
    size[2] = 5;

    typename InputImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);

    typename InputImageType::Pointer image = InputImageType::New();
    image->SetRegions(region);
    image->Allocate();
    image->FillBuffer(0);

    /* Filter */
    using FilterType        = typename itk::CleanLabelSurfaceImageFilter< InputImageType >;
    using FilterPointerType = typename FilterType::Pointer;
    FilterPointerType Filter = FilterType::New();
    Filter->SetInput(image);
    Filter->Update();
    EXPECT_NO_THROW(Filter->Update());

    /* Check */
    itk::ImageRegionIterator< InputImageType > input(Filter->GetOutput(), Filter->GetOutput()->GetLargestPossibleRegion());

    for(input.GoToBegin(); !input.IsAtEnd(); ++input)
    {
        ASSERT_EQ(input.Get(), 0);
    }
}

TEST(itkCleanLabelSurfaceImageFilterUnitTest, RemoveOneSpot) {
    /* Create Image */
    using PixelType  = unsigned int;
    using InputImageType = itk::Image< PixelType, 3 >;

    typename InputImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;

    typename InputImageType::SizeType size;
    size[0] = 5;
    size[1] = 5;
    size[2] = 5;

    typename InputImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);

    typename InputImageType::Pointer image = InputImageType::New();
    image->SetRegions(region);
    image->Allocate();
    image->FillBuffer(0);

    InputImageType::IndexType pixelIndex;
    pixelIndex[0] = 2;
    pixelIndex[1] = 2;
    pixelIndex[2] = 2;
    image->SetPixel(pixelIndex, 1);

    /* Filter */
    using FilterType        = typename itk::CleanLabelSurfaceImageFilter< InputImageType >;
    using FilterPointerType = typename FilterType::Pointer;
    FilterPointerType Filter = FilterType::New();
    Filter->SetInput(image);
    Filter->Update();
    EXPECT_NO_THROW(Filter->Update());

    /* Check */
    itk::ImageRegionIterator< InputImageType > input(Filter->GetOutput(), Filter->GetOutput()->GetLargestPossibleRegion());

    for(input.GoToBegin(); !input.IsAtEnd(); ++input)
    {
        ASSERT_EQ(input.Get(), 0);
    }
}
