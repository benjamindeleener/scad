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
import nibabel as nib
from msct_image import Image
from scad import SCAD
import numpy as np
import scad

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
                print image
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
                manual_segmentation_image.change_orientation()
                centerline_image = Image(image_output)
                centerline_image.change_orientation()

                from msct_types import Coordinate
                # coord_manseg = manual_segmentation_image.getNonZeroCoordinates()
                coord_centerline = centerline_image.getNonZeroCoordinates()

                # check if centerline is in manual segmentation
                result_centerline_in = True
                for coord in coord_centerline:
                    if manual_segmentation_image.data[coord.x, coord.y, coord.z] == 0:
                        result_centerline_in = False
                        print 'failed on slice #' + str(coord.z)
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
            os.chdir('..')

    # t2

    # t2*

    # dmri

    # gre

def validate_scad(folder_input):
    """
    Expecting folder to have the following structure :
    errsm_01:
    - t2
    -- errsm_01.nii.gz or t2.nii.gz
    :param folder_input:
    :return:
    """
    current_folder = os.getcwd()
    os.chdir(folder_input)
    try:
        patients = next(os.walk('.'))[1]
        for i in patients:
            if i != "errsm_01" and i !="errsm_02":
                directory = i + "/t2"
                os.chdir(directory)
                try:
                    if os.path.isfile(i+"_t2.nii.gz"):
                        raw_image = Image(i+"_t2.nii.gz")
                    elif os.path.isfile("t2.nii.gz"):
                        raw_image = Image("t2.nii.gz")
                    else:
                        raise Exception("t2.nii.gz or "+i+"_t2.nii.gz file is not found")

                    raw_orientation = raw_image.change_orientation()
                    SCAD(raw_image, contrast="t2", rm_tmp_file=1, verbose=1).test_debug()

                    manual_seg = Image(i+"_t2_manual_segmentation.nii.gz")
                    manual_orientation = manual_seg.change_orientation()

                    from scipy.ndimage.measurements import center_of_mass
                    # find COM
                    iterator = range(manual_seg.data.shape[2])
                    com_x = [0 for ix in iterator]
                    com_y = [0 for iy in iterator]

                    for iz in iterator:
                        com_x[iz], com_y[iz] = center_of_mass(manual_seg.data[:, :, iz])
                    #raw_image.change_orientation(raw_orientation)
                    #manual_seg.change_orientation(manual_orientation)

                    centerline_scad = Image(i+"_t2_centerline.nii.gz")
                    os.remove(i+"_t2_centerline.nii.gz")

                    centerline_scad.change_orientation()
                    distance = []
                    for iz in range(centerline_scad.data.shape[2]):
                        ind1 = np.argmax(centerline_scad.data[:, :, iz])
                        X,Y = scad.ind2sub(centerline_scad.data[:, :, i].shape,ind1)
                        com_phys = centerline_scad.transfo_pix2phys([[com_x[iz], com_y[iz], iz]])
                        scad_phys = centerline_scad.transfo_pix2phys([[X, Y, iz]])
                        distance_magnitude = np.linalg.norm(com_phys-scad_phys)
                        distance.append(distance_magnitude)



                    os.chdir(folder_input)

                except Exception, e:
                    print e.message
                pass
    except Exception, e:
        print e.message



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

    ### Start of not good code
    if "-scad" in script_arguments:
        folder = script_arguments[script_arguments.index("-i") + 1]
        if folder != "" or folder is not None:
            validate_scad(folder)
    # elif len(script_arguments) > 3:
    #     print 'ERROR: this script only accepts three arguments: the name of your algorithm, if it is a python script or' \
    #           'not and the verbose option.'
    #     usage()
    #
    # algorithm = script_arguments[0]
    # verbose = True
    # ispython = False
    # if len(script_arguments) >= 2:
    #     if 'verbose' in script_arguments[1:]:
    #         verbose = False
    #     if 'ispython' in script_arguments[1:]:
    #         ispython = True
    #
    # scadMRValidation(algorithm=algorithm, isPython=ispython, verbose=verbose)
