#!/usr/bin/env python
#########################################################################################
#
# Validation script for SCAD (Spinal Cord Automatic Detection)
#
# Brainhack MTL 2015: Algorithms for automatic spinal cord detection on MR images
#
# This repository is intented to develop and test new algorithms for automatically detect the spinal cord on various
# contrasts of MR volumes.
# The developed algorithms must satisfy the following criteria:
# - they can be coded in Python or C++
# - they must read a nifti image as input image (.nii or .nii.gz): "-i" (input file name) option
# - they have to output a binary image with the same format and orientation as input image, containing the location
#   or the centerline of the spinal cord: "-o" (output file name) option
# - they have to be **fast**
#
# To validate a new algorithm, it must go through the validation pipeline using the following command:
#
# scad_validation.py "algo_name"
#
# The validation pipeline tests your algorithm throughout a testing dataset, containing many images of the spinal cord
# with various contrasts and fields of view, along with their manual segmentation.
# It tests several criteria:
# 1. if your detection is inside the spinal cord
# 2. if your detection is near the spinal cord centerline (at least near the manual centerline)
# 3. if the length of the centerline your algorithm extracted correspond with the manual centerline
#
# ---------------------------------------------------------------------------------------
# Copyright (c) 2015 Polytechnique Montreal <www.neuro.polymtl.ca>
# Authors: Benjamin De Leener
# Modified: 2015-07-22
#
# About the license: see the file LICENSE
#########################################################################################

import sys
import os

def scadMRValidation(algorithm, isPython=False, verbose=True):
    if not isinstance(algorithm, str) or not algorithm:
        print 'ERROR: You must provide the name of your algorithm as a string.'
        usage()

    import time
    import sct_utils as sct

    # creating a new folder with the experiment
    path_experiment = 'scad-experiment.'+algorithm+'.'+time.strftime("%y%m%d%H%M%S")
    #status, output = sct.run('mkdir '+path_experiment, verbose)

    # copying images from "data" folder into experiment folder
    sct.copyDirectory('data', path_experiment)

    # Starting validation
    os.chdir(path_experiment)
    # t1
    os.chdir('t1/')
    for subject_dir in os.listdir('./'):
        if os.path.isdir(subject_dir):
            os.chdir(subject_dir)

            # creating list of images and corresponding manual segmentation
            list_images = dict()
            for file_name in os.listdir('./'):
                if not 'manual_segmentation' in file_name:
                    for file_name_corr in os.listdir('./'):
                        if 'manual_segmentation' in file_name_corr and sct.extract_fname(file_name)[1] in file_name_corr:
                            list_images[file_name] = file_name_corr

            # running the proposed algorithm on images in the folder and analyzing the results
            for image, image_manual_seg in list_images.items():
                path_in, file_in, ext_in = sct.extract_fname(image)
                image_output = file_in+'_centerline'+ext_in
                if ispython:
                    try:
                        eval(algorithm+'('+image+', t1, verbose='+str(verbose)+')')
                    except Exception as e:
                        print 'Error during spinal cord detection on line {}:'.format(sys.exc_info()[-1].tb_lineno)
                        print 'Subject: t1/'+subject_dir+'/'+image
                        print e
                        sys.exit(2)
                else:
                    cmd = algorithm+' -i '+image+' -t t1'
                    if verbose:
                        cmd += ' -v'
                    status, output = sct.run(cmd, verbose=verbose)
                    if status != 0:
                        print 'Error during spinal cord detection on Subject: t1/'+subject_dir+'/'+image
                        print output
                        sys.exit(2)

                # analyzing the resulting centerline
                from msct_image import Image
                manual_segmentation_image = Image(image_manual_seg)
                centerline_image = Image(image_output)

                from msct_types import Coordinate
                # coord_manseg = manual_segmentation_image.getNonZeroCoordinates()
                coord_centerline = centerline_image.getNonZeroCoordinates()

                # check if centerline is in manual segmentation
                result_centerline_in = True
                for coord in coord_centerline:
                    if manual_segmentation_image.data[coord.x, coord.y, coord.z] == 0:
                        result_centerline_in = False
                        break
                if result_centerline_in:
                    print 'OK: Centerline is inside manual segmentation.'
                else:
                    print 'FAIL: Centerline is outside manual segmentation.'

                # check the length of centerline compared to manual segmentation
                # import sct_process_segmentation as sct_seg
                # length_manseg = sct_seg.compute_length(image_manual_seg)
                # length_centerline = sct_seg.compute_length(image_output)
                # if length_manseg*0.9 <= length_centerline <= length_manseg*1.1:
                #     print 'OK: Length of centerline correspond to length of manual segmentation.'
                # else:
                #     print 'FAIL: Length of centerline does not correspond to length of manual segmentation.'

    # t2

    # t2*

    # dmri

    # gre


def usage():
    print """
    """ + os.path.basename(__file__) + """
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Brainhack MTL 2015

    DESCRIPTION
      Validation script for SCAD (Spinal Cord Automatic Detection)

    USAGE
      """ + os.path.basename(__file__) + """ <algorithm_name>

    MANDATORY ARGUMENTS
      <algorithm_name>  name of the script you want to validate. The script must have -i, -o and -v options enabled.

    OPTIONAL ARGUMENTS
      ispython          Switch to python validation. It means that the algorithm will be called as a python method.
      verbose           Disable display. Default: display on.
      -h                help. Show this message
    """
    sys.exit(1)

# START PROGRAM
# ==========================================================================================
if __name__ == "__main__":
    # reading the name of algorithm from arguments
    script_arguments = sys.argv[1:]
    if "-h" in script_arguments:
        usage()
    elif len(script_arguments) > 3:
        print 'ERROR: this script only accepts three arguments: the name of your algorithm, if it is a python script or' \
              'not and the verbose option.'
        usage()

    algorithm = script_arguments[0]
    verbose = True
    ispython = False
    if len(script_arguments) >= 2:
        if 'verbose' in script_arguments[1:]:
            verbose = False
        if 'ispython' in script_arguments[1:]:
            ispython = True

    scadMRValidation(algorithm=algorithm, isPython=ispython, verbose=verbose)
