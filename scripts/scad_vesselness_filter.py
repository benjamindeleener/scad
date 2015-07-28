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
import sct_utils as sct
import sys
import time
import os
import numpy as np
# import dipy

class Param:
    def __init__(self):
        self.debug = 0
        self.verbose = 1


class Vesselness:
    def __init__(self, fname):
        self.image = Image(fname)

        self.pretreated_im = self.pretreat()

        self.result_im = self.image.copy()
        self.result_im.file_name = self.image.file_name + '_vessel_mask'
        self.result_im.path = './'

        self.hess = hessian(self.pretreated_im.data)
        self.vessel_mask = self.compute_vessel_mask()

        self.result_im.data = self.vessel_mask
        self.result_im.save()

    def pretreat(self):
        import scipy.ndimage.filters as scp_filters
        '''
        status, orientation = sct.run('sct_orientation.py -i ' + self.image.path + self.image.file_name + self.image.ext)
        orientation = orientation[4:7]
        if orientation != 'RPI':
            sct.run('sct_orientation.py -i ' + self.image.path + self.image.file_name + self.image.ext + ' -s RPI')
        fname = self.image.file_name[:-7] + '_RPI.nii.gz'
        '''
        fname = self.image.file_name

        pretreated = self.image.copy()
        pretreated.file_name = fname

        nx, ny, nz, nt, px, py, pz, pt = sct.get_dimension(pretreated.file_name)
        sc_size = 3  # in mm
        sc_npix = ((sc_size/px) + (sc_size/py))/2.0
        pretreated.data = scp_filters.gaussian_filter(pretreated.data, sigma=sc_npix)
        pretreated.file_name = pretreated.file_name + '_gaussian'
        pretreated.save()
        return pretreated


    def compute_vessel_mask(self):
        print 'COMPUTE VESSEL MASK'
        vessel_mask = np.zeros(self.pretreated_im.data.shape)
        print 'Image shape : ', str(self.pretreated_im.data.shape)
        for x in range(self.pretreated_im.data.shape[0]):
            print '----------> x =', x
            now_x = time.time()
            for y in range(self.pretreated_im.data.shape[1]):
                for z in range(self.pretreated_im.data.shape[2]):
                    # H = self.hess[x][y][z]
                    H = self.get_hessian(x, y, z)
                    vals, vects = np.linalg.eig(H)
                    l1, l2, l3 = sort_eigen_vals(vals)
                    vessel_mask[x][y][z] = line_filter(l1, l2, l3)
            print 'x ' + str(x) + ' --> in ' + str(time.time() - now_x) + ' sec'
        return vessel_mask

    def get_hessian(self, x, y, z):
        Hxx = self.hess[0][0][x][y][z]
        Hxy = self.hess[0][1][x][y][z]
        Hxz = self.hess[0][2][x][y][z]
        Hyx = self.hess[1][0][x][y][z]
        Hyy = self.hess[1][1][x][y][z]
        Hyz = self.hess[1][2][x][y][z]
        Hzx = self.hess[2][0][x][y][z]
        Hzy = self.hess[2][1][x][y][z]
        Hzz = self.hess[2][2][x][y][z]
        H = [[Hxx, Hxy, Hxz],
             [Hyx, Hyy, Hyz],
             [Hzx, Hzy, Hzz]]
        return H


def hessian(im_data):
    print 'COMPUTE HESSIAN'
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
                Hxx = hess[0][0][x][y][z]
                Hxy = hess[0][1][x][y][z]
                Hxz = hess[0][2][x][y][z]
                Hyx = hess[1][0][x][y][z]
                Hyy = hess[1][1][x][y][z]
                Hyz = hess[1][2][x][y][z]
                Hzx = hess[2][0][x][y][z]
                Hzy = hess[2][1][x][y][z]
                Hzz = hess[2][2][x][y][z]
                im_hess[x][y][z] = [[Hxx, Hxy, Hxz],
                                    [Hyx, Hyy, Hyz],
                                    [Hzx, Hzy, Hzz]]
    '''
    return hess  # im_hess


def line_filter(lambda1, lambda2, lambda3, alpha1=0.1,  alpha2=1.0, alpha3=5.0):
    """
    detect a bright tubular structure
    eq 1 from http://www.spl.harvard.edu/archive/spl-pre2007/pages/papers/yoshi/node3.html#SECTION00021000000000000000
    eq 2 from Generalizing vesselness with respect to dimensionality and shape, Luca Antiga - 2007
    :param lambda1:
    :param lambda2:
    :param lambda3:
    :return:
    """
    from math import exp, sqrt
    '''
    # eq 1
    lambda_c = min(-lambda2, -lambda3)

    if lambda_c == 0:
        return 0
    elif lambda1 > 0:
        return lambda_c * exp(-lambda1**2/2*(alpha2*lambda_c)**2)
    elif lambda1 <= 0:
        return lambda_c * exp(-lambda1**2/2*(alpha1*lambda_c)**2)
    '''

    # eq 2
    if lambda2 < 0 and lambda3 < 0:
        RA = abs(lambda3)/abs(lambda3)
        RB = abs(lambda1)/abs(lambda2*lambda3)
        S = sqrt(lambda1**2 + lambda2**2 + lambda3**2)
        return (1-exp(-RA**2/2*alpha1**2))*exp(-RB**2/2*alpha2**2)*(1-exp(-S**2/2*alpha3**2))
    else:
        return 0


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
    tmp_dir = 'tmp.' + time.strftime("%y%m%d%H%M%S")
    sct.run('mkdir ' + tmp_dir)
    im_name = 'input_image.nii.gz'
    sct.run('cp ' + input_filename + ' ./' + tmp_dir + '/' + im_name)

    os.chdir(tmp_dir)
    vessel_filter = Vesselness(im_name)
    sct.run('cp ' + vessel_filter.result_im.file_name + vessel_filter.result_im.ext + '  ../' + sct.extract_fname(input_filename)[1] + '_vessel_filtered.nii.gz')
    os.chdir('..')