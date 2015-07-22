# Spinal Cord Automatic Detection
Brainhack MTL 2015: Algorithms for automatic spinal cord detection on MR images

This repository is intented to develop and test new algorithms for automatically detect the spinal cord on various contrasts of MR volumes.
The developed algorithms must satisfy the following criteria:
- they can be coded in Python or C++
- they must read a nifti image as input image (.nii or .nii.gz): "-i" (input file name) option
- they have to output a binary image with the same format and orientation as input image, containing the location or the centerline of the spinal cord: "-o" (output file name) option
- they have to be **fast**

To validate a new algorithm, it must go through the validation pipeline using the following command:
```
scad_validation.py "algo_name"
```

The validation pipeline tests your algorithm throughout a testing dataset, containing many images of the spinal cord with various contrasts and fields of view, along with their manual segmentation.
It tests several criteria:

1. if your detection is inside the spinal cord
2. if your detection is near the spinal cord centerline (at least near the manual centerline)
3. if the length of the centerline your algorithm extracted correspond with the manual centerline


## License

The MIT License (MIT)

Copyright (c) 2015 Benjamin De Leener and contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
