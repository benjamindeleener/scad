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
import os
import sct_utils as sct


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
        # create tmp and copy input
        path_tmp = sct.tmp_create()
        sct.tmp_copy_nifti(self.input_image.absolutepath, path_tmp, 'raw.nii')
        os.chdir(path_tmp)

        import glob
        matlab_run=glob.glob('/Applications/MATLAB*')[0]+'/bin/matlab -nodesktop -nosplash -r'

        # SymmetryDetector
        sct.run(matlab_run + ' "scad_symetry(\'raw.nii\'); exit"')
        # minimal path
        #sct.run(matlab_run + ' "scad_minimalpath_homo(\'raw.nii\'); exit"')
        # vesselness filter
        sct.run('sct_vesselness -i raw.nii -t ' + self._contrast)
        # minimal path on vesselness filter
        sct.run(matlab_run + ' "scad_minimalpath(\'imageVesselNessFilter.nii.gz\'); exit"')
        # multiply results
        sct.run('fslmaths imageVesselNessFilter_minimalpath.nii -mul raw_symetry.nii.gz spine_detection.nii.gz')

        # get centerline
        sct.run(matlab_run + ' "scad_proba2centerline(\'spine_detection.nii.gz\'); exit"')

        # copy back centerline
        os.chdir('../')
        sct.tmp_copy_nifti(path_tmp + 'centerline.nii.gz',self.input_image.path,self.input_image.file_name+'_centerline'+self.input_image.ext)

        # delete tmp
        sct.run('rm -rf ' + path_tmp)
        # symDetector = SymmetryDetector(self.contrast, self.verbose)
        # symDetector.execute()
        # SymmetryDetector()



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


if __name__ == "__main__":
    parser = Script.get_parser()

    arguments = parser.parse(sys.argv[1:])

    input_image = Image(arguments["-i"])
    contrast_type = arguments["-t"]

    scad = SCAD(input_image)
    scad.contrast = contrast_type
    scad.execute()
