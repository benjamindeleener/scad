#!/usr/bin/env python
#########################################################################################
#
# Spinal Cord Automatic Detection
#
# ---------------------------------------------------------------------------------------
# Copyright (c) 2015 Polytechnique Montreal <www.neuro.polymtl.ca>
# Authors: Benjamin De Leener
# Modified: 2015-07-27
#
# About the license: see the file LICENSE
#########################################################################################

import sys
from msct_base_classes import BaseScript, Algorithm
from msct_parser import Parser
from msct_image import Image


class Script(BaseScript):
    def __init__(self):
        super(Script, self).__init__()

    @staticmethod
    def get_parser():
        # Initialize the parser
        parser = Parser(__file__)
        parser.usage.set_description('''This program automatically detect the spinal cord in a MR image and output a centerline of the spinal cord.''')
        parser.add_option(name="-i",
                          type_value="file",
                          description="input image.",
                          mandatory=True,
                          example="t2.nii.gz")
        parser.add_option(name="-t",
                          type_value="multiple_choice",
                          description="type of image contrast, t2: cord dark / CSF bright ; t1: cord bright / CSF dark",
                          mandatory=True,
                          example=['t1', 't2'])
        parser.usage.addSection("General options")
        parser.add_option(name="-v",
                          type_value="multiple_choice",
                          description="1: display on, 0: display off (default)",
                          mandatory=False,
                          example=["0", "1"],
                          default_value="1")
        parser.add_option(name="-h",
                          type_value=None,
                          description="display this help",
                          mandatory=False)

        return parser



class SCAD(Algorithm):
    def __init__(self, input_image, contrast=None, verbose=0):
        super(SCAD, self).__init__(input_image)
        self._contrast = contrast
        self._verbose = verbose

    @property
    def contrast(self):
        return self._contrast

    @contrast.setter
    def contrast(self, value):
        if value in ['t1', 't2']:
            self._contrast = value
        else:
            raise Exception('ERROR: contrast value must be t1 or t2')

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, value):
        if value in [0, 1]:
            self._verbose = value
        else:
            raise Exception('ERROR: verbose value must be an integer and equal to 0 or 1')

    def execute(self):
        print 'Execution of the SCAD algorithm'
        symDetector = SymmetryDetector(self.contrast, self.verbose)
        symDetector.execute()


class SymmetryDetector(Algorithm):
    def __init__(self, input_image, contrast=None, verbose=0):
        super(SymmetryDetector, self).__init__(input_image)
        self._contrast = contrast
        self._verbose = verbose

    @property
    def contrast(self):
        return self._contrast

    @contrast.setter
    def contrast(self, value):
        if value in ['t1', 't2']:
            self._contrast = value
        else:
            raise Exception('ERROR: contrast value must be t1 or t2')

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, value):
        if value in [0, 1]:
            self._verbose = value
        else:
            raise Exception('ERROR: verbose value must be an integer and equal to 0 or 1')

    def execute(self):
        print 'Execution of the symmetry detection algorithm'

        # for each slice, cut the slice in two at the middle point in the left-right direction




        map < double, int > mutualInformation;
        int
        startSlice = desiredSizeInitial[2] / 4, endSlice = desiredSizeInitial[2] / 4 * 3;
        if (desiredSizeInitial[2] < cropSize * 2) {
        startSlice = cropSize / 2;
        endSlice = desiredSizeInitial[2]-cropSize / 2;
        }

    // Check for non - null intensity in the image.If null, mutual information cannot be computed...
    ImageType::IndexType
    desiredStart_i;


ImageType::SizeType
desiredSize_i = desiredSizeInitial;
desiredStart_i[0] = 0;
desiredStart_i[1] = (int)
init_slice;
desiredStart_i[2] = 0;
desiredSize_i[1] = 0;
desiredSize_i[2] = desiredSizeInitial[2];
ImageType::RegionType
desiredRegionImage(desiredStart_i, desiredSize_i);
typedef
itk::ExtractImageFilter < ImageType, ImageType2D > Crop2DFilterType;
Crop2DFilterType::Pointer
cropFilter = Crop2DFilterType::New();
cropFilter->SetExtractionRegion(desiredRegionImage);
cropFilter->SetInput(inputImage_);
# if ITK_VERSION_MAJOR >= 4
cropFilter->SetDirectionCollapseToIdentity(); // This is required.
# endif
try {
cropFilter->Update();
} catch( itk::
ExceptionObject & e ) {
    std::cerr << "Exception caught while updating cropFilter " << std::endl;
std::cerr << e << std::endl;
}
ImageType2D::Pointer
image = cropFilter->GetOutput();
typedef
itk::MinimumMaximumImageCalculator < ImageType2D > MinMaxCalculatorType;
MinMaxCalculatorType::Pointer
minMaxCalculator = MinMaxCalculatorType::New();
minMaxCalculator->SetImage(image);
minMaxCalculator->ComputeMaximum();
minMaxCalculator->ComputeMinimum();
ImageType2D::PixelType
maxIm = minMaxCalculator->GetMaximum(), minIm = minMaxCalculator->GetMinimum();
if (maxIm == minIm) {
cerr << "ERROR: The axial slice where the symmetry will be detected (slice " << init_slice << ") is full of constant value (" << maxIm << "). You can change it using -init parameter." << endl;
return -1;
}

for (int i=startSlice; i < endSlice; i++)
{
float
startCrop = i, size;
if (startCrop < desiredSizeInitial[2] / 2 & & startCrop <= bandWidth_ + 1) size = startCrop-1;
else if (startCrop >= desiredSizeInitial[2] / 2 & & startCrop >= desiredSizeInitial[2]-bandWidth_-1) size = desiredSizeInitial[2]-startCrop-1;
else size = bandWidth_;
ImageType::IndexType
desiredStart;
ImageType::SizeType
desiredSize = desiredSizeInitial;
desiredStart[0] = 0;
desiredStart[1] = (int)
init_slice;
desiredStart[2] = startCrop;
desiredSize[1] = 0;
desiredSize[2] = size;

// Right
Image
ImageType::RegionType
desiredRegionImageRight(desiredStart, desiredSize);
typedef
itk::ExtractImageFilter < ImageType, ImageType2D > Crop2DFilterType;
Crop2DFilterType::Pointer
cropFilterRight = Crop2DFilterType::New();
cropFilterRight->SetExtractionRegion(desiredRegionImageRight);
cropFilterRight->SetInput(inputImage_);
# if ITK_VERSION_MAJOR >= 4
cropFilterRight->SetDirectionCollapseToIdentity(); // This is required.
# endif
try {
cropFilterRight->Update();
} catch( itk::
ExceptionObject & e ) {
    std::cerr << "Exception caught while updating cropFilter " << std::endl;
std::cerr << e << std::endl;
}
ImageType2D::Pointer
imageRight = cropFilterRight->GetOutput();

// Left
Image
desiredStart[2] = startCrop - size;
if (desiredStart[2] < 0) desiredStart[2] = 0;
ImageType::RegionType
desiredRegionImageLeft(desiredStart, desiredSize);
Crop2DFilterType::Pointer
cropFilterLeft = Crop2DFilterType::New();
cropFilterLeft->SetExtractionRegion(desiredRegionImageLeft);
cropFilterLeft->SetInput(inputImage_);
# if ITK_VERSION_MAJOR >= 4

cropFilterLeft->SetDirectionCollapseToIdentity(); // This is required.
# endif
try {
cropFilterLeft->Update();
} catch( itk::
ExceptionObject & e ) {
    std::cerr << "Exception caught while updating cropFilter " << std::endl;
std::cerr << e << std::endl;
}
ImageType2D::Pointer
imageLeft = cropFilterLeft->GetOutput();

ImageType2D::IndexType
desIndex;
desIndex.Fill(0);
ImageType2D::SizeType
desSize;
desSize[0] = desiredSize[0];
desSize[1] = desiredSize[2];
ImageType2D::RegionType
desired2DRegionImageRight(desIndex, desSize);
ImageType2D::RegionType
desired2DRegionImageLeft(desIndex, desSize);
imageRight->SetLargestPossibleRegion(desired2DRegionImageRight);
imageRight->SetRequestedRegion(desired2DRegionImageRight);
imageRight->SetRegions(desired2DRegionImageRight);
imageLeft->SetLargestPossibleRegion(desired2DRegionImageLeft);
imageLeft->SetRequestedRegion(desired2DRegionImageLeft);
imageLeft->SetRegions(desired2DRegionImageLeft);

itk::FixedArray < bool, 2 > flipAxes;
flipAxes[0] = false;
flipAxes[1] = true;
typedef
itk::FlipImageFilter < ImageType2D > FlipImageFilterType;
FlipImageFilterType::Pointer
flipFilter = FlipImageFilterType::New();
flipFilter->SetInput(imageRight);
flipFilter->SetFlipAxes(flipAxes);
flipFilter->Update();
imageRight = flipFilter->GetOutput();

ImageType2D::PointType
origin = imageLeft->GetOrigin();
imageRight->SetOrigin(origin);

// Better
value is minimum
typedef
itk::MattesMutualInformationImageToImageMetric < ImageType2D, ImageType2D > MattesMutualInformationFilter;
MattesMutualInformationFilter::Pointer
correlationFilter = MattesMutualInformationFilter::New();
typedef
itk::IdentityTransform < double, 2 > IdentityTransform;
IdentityTransform::Pointer
transform = IdentityTransform::New();
correlationFilter->SetTransform(transform);
typedef
itk::NearestNeighborInterpolateImageFunction < ImageType2D, double > InterpolatorType;
InterpolatorType::Pointer
interpolator = InterpolatorType::New();
interpolator->SetInputImage(imageRight);
correlationFilter->SetInterpolator(interpolator);
correlationFilter->SetFixedImage(imageLeft);
correlationFilter->SetMovingImage(imageRight);
correlationFilter->SetFixedImageRegion(imageLeft->GetLargestPossibleRegion());
correlationFilter->UseAllPixelsOn();
correlationFilter->Initialize();
MattesMutualInformationFilter::TransformParametersType
id(3);

id[0] = 0;
id[1] = 0;
id[2] = 0;
double
value = 0.0;
try {
value = correlationFilter->GetValue( id );
} catch( itk::
ExceptionObject & e ) {
    std::cerr << "Exception caught while getting value " << std::endl;
std::cerr << e << std::endl;
}
mutualInformation[value] = startCrop;
}
// cout << "Cropping around slice = " << mutualInformation.begin()->second << endl;
middleSlice_ = mutualInformation.begin()->second;
return middleSlice_;

if __name__ == "__main__":
    parser = Script.get_parser()

    arguments = parser.parse(sys.argv[1:])

    input_image = Image(arguments["-i"])
    contrast_type = arguments["-t"]

    scad = SCAD(input_image)
    scad.contrast = contrast_type
    scad.execute()
