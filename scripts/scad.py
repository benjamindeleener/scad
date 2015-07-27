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
    def __init__(self):
        super(SCAD, self).__init__()
        self._contrast = None

    @property
    def contrast(self):
        return self._contrast

    @contrast.setter
    def contrast(self, value):
        if value in ['t1', 't2']:
            self._contrast = value
        else:
            raise Exception('ERROR: contrast value must be t1 or t2')

    def execute(self):
        print 'Execution of the algorithm'



if __name__ == "__main__":
    parser = Script.get_parser()

    arguments = parser.parse(sys.argv[1:])

    input_image = Image(arguments["-i"])
    contrast_type = arguments["-t"]

    scad = SCAD()
    scad.input_image = input_image
    scad.contrast = contrast_type
    scad.execute()