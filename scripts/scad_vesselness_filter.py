#!/usr/bin/env python
#########################################################################################
#
# Vesselness filter
# hessian computation inpired from: https://code.google.com/p/pycvf/source/browse/pycvf/trunk/pycvf/nodes/image/hessian.py?spec=svn22&r=22
#
# ---------------------------------------------------------------------------------------
# Copyright (c) 2015 Polytechnique Montreal <www.neuro.polymtl.ca>
# Author: Sara Dupont
# Modified: 2015-07-27
#
# About the license: see the file LICENSE.TXT
#########################################################################################
from msct_parser import Parser
from msct_image import Image
import sys
import numpy as np
# import dipy

class Param:
    def __init__(self):
        self.debug = 0
        self.verbose = 1


class Vesselness:
    def __init__(self, fname):
        self.image = Image(fname)
        self.hess = hessian(self.image.data)
        self.vessel_mask = self.compute_vessel_mask()

    def compute_vessel_mask(self):
        vessel_mask = np.zeros(self.image.data.shape)
        for x in range(self.image.data.shape[0]):
            for y in range(self.image.data.shape[1]):
                for z in range(self.image.data.shape[2]):

                    Hxx = self.hess[0][0].tolist()[x][y][z]
                    Hxy = self.hess[0][1].tolist()[x][y][z]
                    Hxz = self.hess[0][2].tolist()[x][y][z]
                    Hyx = self.hess[1][0].tolist()[x][y][z]
                    Hyy = self.hess[1][1].tolist()[x][y][z]
                    Hyz = self.hess[1][2].tolist()[x][y][z]
                    Hzx = self.hess[2][0].tolist()[x][y][z]
                    Hzy = self.hess[2][1].tolist()[x][y][z]
                    Hzz = self.hess[2][2].tolist()[x][y][z]
                    H = [[Hxx, Hxy, Hxz],
                        [Hyx, Hyy, Hyz],
                        [Hzx, Hzy, Hzz]]

                    # H = self.hess[x][y][z]
                    vals, vects = np.linalg.eig(H)
                    l1, l2, l3 = sort_eigen_vals(vals)
                    vessel_mask[x][y][z] = line_filter(l1, l2, l3)
        return vessel_mask


def hessian(im_data):
    grad = np.gradient(im_data)  # grad[0] = X gradient, grad[1] = Y gradient, grad[2] = Z gradient
    hess = []
    for grad_i in grad:
        hess.append(np.gradient(grad_i))
    '''
    import dipy.align.metrics as dpm
    return dpm.gradient(dpm.gradient(im_data))
    '''
    # hess[0][0] = XX Hessian, ...

    '''
    im_hess = np.zeros(im_data.shape).tolist()
    for x in range(im_data.shape[0]):
        for y in range(im_data.shape[1]):
            for z in range(im_data.shape[2]):
                Hxx = hess[0][0].tolist()[x][y][z]
                Hxy = hess[0][1].tolist()[x][y][z]
                Hxz = hess[0][2].tolist()[x][y][z]
                Hyx = hess[1][0].tolist()[x][y][z]
                Hyy = hess[1][1].tolist()[x][y][z]
                Hyz = hess[1][2].tolist()[x][y][z]
                Hzx = hess[2][0].tolist()[x][y][z]
                Hzy = hess[2][1].tolist()[x][y][z]
                Hzz = hess[2][2].tolist()[x][y][z]
                im_hess[x][y][z] = [[Hxx, Hxy, Hxz],
                                    [Hyx, Hyy, Hyz],
                                    [Hzx, Hzy, Hzz]]
    '''
    return hess  # im_hess


def line_filter(lambda1, lambda2, lambda3, alpha1=0.5, alpha2=2):
    """
    detect a bright tubular structure
    from http://www.spl.harvard.edu/archive/spl-pre2007/pages/papers/yoshi/node3.html#SECTION00021000000000000000
    :param lambda1:
    :param lambda2:
    :param lambda3:
    :return:
    """
    from math import exp
    lambda_c = min(-lambda2, -lambda3)

    if lambda_c == 0:
        return 0
    elif lambda1 > 0:
        return lambda_c * exp(-lambda1**2/2*(alpha2*lambda_c)**2)
    elif lambda1 <= 0:
        return lambda_c * exp(-lambda1**2/2*(alpha1*lambda_c)**2)


def sort_eigen_vals(l):
    l_abs = list(abs(np.asarray(l)))
    d = zip(l_abs, range(len(l_abs)))
    sorted_l = []
    for v, i in d:
        sorted_l.append(l[i])
    return sorted_l

'''
def hessian(im):
    im_dim = len(im.shape)
    if im_dim == 2:
        res = np.zeros((im_dim, im_dim), dtype=object)
        for i in range(im_dim):
            for j in range(im_dim):
                res[i, j] = np.diff(np.diff(im, axis=j), axis=i)
    elif im_dim == 3:
        res = np.zeros((im_dim, im_dim, im_dim), dtype=object)
        for i in range(im_dim):
            for j in range(im_dim):
                for k in range(im_dim):
                    res[i, j, k] = np.diff(np.diff(im, axis=j), axis=i) #TODO: adapt for 3D
    return res
'''



#=======================================================================================================================
# Start program
#=======================================================================================================================
if __name__ == "__main__":
    # initialize parameters
    param = Param()
    param_default = Param()

    # Initialize the parser
    parser = Parser(__file__)
    parser.usage.set_description('Vesselness filter - in developement')
    parser.add_option(name="-i",
                      type_value="file",
                      description="Image to filter.",
                      mandatory=True,
                      example="t2.nii.gz")
    parser.add_option(name="-v",
                      type_value="multiple_choice",
                      description="verbose. Default=" + str(param_default.verbose),
                      mandatory=False,
                      example=['0', '1', '2'])
    arguments = parser.parse(sys.argv[1:])

    input_filename = arguments["-i"]

    if "-v" in arguments:
        param.verbose = arguments["-v"]

    vessel_filter = Vesselness(input_filename)
    Image(param=vessel_filter.vessel_mask, absolutepath='test_vessel_mask.nii.gz').save()