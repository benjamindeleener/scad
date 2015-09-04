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
from scipy.ndimage import interpolation
from scipy.stats._continuous_distns import invgamma_gen
from msct_base_classes import BaseScript, Algorithm
from msct_parser import Parser
from msct_image import Image
import os
import sct_utils as sct
import numpy as np
import nibabel as nb
from sct_straighten_spinalcord import smooth_centerline
from skimage.morphology import watershed


def treshold_min_path(data, threshold):
    result = data
    result[result >= threshold] = 1
    result[result < threshold] = 0

    return result


def watershed_segmentation(img_data):
    pass

def get_line_cross_corr():
    pass


def block_symmetry(img, section_nb):
    """
    This function splits an image in section_nb number of blocks along the A-P axis and calculates the l-r symmetry for each of them. Then, concatenantes them to provide the symmetry of the image
    :param img: the image
    :param section_nb: number of sections that the img will be divided into
    :return: Concatenated sections, each having its own symmetry
    """
    raw_orientation = img.change_orientation()
    dim = img.data.shape

    section_length = dim[1]/section_nb

    result = np.zeros(dim)

    for i in range(0, section_nb):
        if (i+1)*section_length > dim[1]:
            y_length = (i+1)*section_length - ((i+1)*section_length - dim[1])
            result[:, i*section_length:i*section_length + y_length, :] = symmetry_detector_right_left(img.data[:, i*section_length:i*section_length + y_length, :], (img.data[:, i*section_length:i*section_length + y_length, :]).shape, cropped_xy=1)
        sym = symmetry_detector_right_left(img.data[:, i*section_length:(i+1)*section_length, :], (img.data[:, i*section_length:(i+1)*section_length, :]).shape, cropped_xy=1)
        result[:, i*section_length:(i+1)*section_length, :] = sym

    return result


def mixed_symmetry():
    pass


def equalize_array_histogram(array):
    array_min = np.amin(array)
    array -= array_min
    array_max = np.amax(array)
    array /= array_max

    return array


def get_minimum_path(data, smooth_factor=np.sqrt(2), invert=0, verbose=1):
    [m, n, p] = data.shape
    max_value = np.amax(data)
    if invert:
        data=max_value-data
    J1 = np.ones([m, n, p])*np.inf
    J2 = J1
    J1[:, :, 0] = 0
    for row in range(0, p-1):
        pJ = J1[:, :, row-1]
        cP = np.squeeze(data[1:-2, 1:-2, row])
        # cP = np.squeeze(data[:, :, row])
        VI = np.dstack((cP*smooth_factor, cP*smooth_factor, cP, cP*smooth_factor, cP*smooth_factor))

        Jq = np.dstack((pJ[0:-3, 1:-2], pJ[2:-1, 1:-2], pJ[1:-2, 1:-2], pJ[1:-2, 0:-3], pJ[1:-2, 2:-1]))
        J1[1:-2, 1:-2, row] = np.min(Jq+VI, 2)
        pass

    J2[:, :, p-1] = 0
    for row in range(p-2, 1, -1):
        pJ = J2[:, :, row+1]
        cP = np.squeeze(data[1:-2, 1:-2, row])
        # cP = np.squeeze(data[:, :, row])
        VI = np.dstack((cP*smooth_factor, cP*smooth_factor, cP, cP*smooth_factor, cP*smooth_factor))

        Jq = np.dstack((pJ[0:-3, 1:-2], pJ[2:-1, 1:-2], pJ[1:-2,1:-2], pJ[1:-2, 0:-3], pJ[1:-2, 2:-1]))
        J2[1:-2, 1:-2, row] = np.min(Jq+VI, 2)
        pass

    result = J1+J2
    if invert:
        percent = np.percentile(result, 50)
        result[result>50] = percent

        result_min = np.amin(result)
        result_max = np.amax(result)
        result = np.divide(np.subtract(result, result_min), result_max)
        result_max = np.amax(result)


    for i in range(1, p-2):
        result[:, :, i] = 1 - result[:,:,i]

    return result


def get_minimum_path_nii(fname):
    from msct_image import Image
    data=Image(fname)
    vesselness_data = data.data
    raw_orient=data.change_orientation()
    data.data=get_minimum_path(data.data, invert=1)
    data.change_orientation(raw_orient)
    data.file_name += '_minimalpath'
    data.save()


def remove_outliers(img_centerline_unsmoothed, std=None,mu=None ):
    centerline_unsmoothed = img_centerline_unsmoothed.data
    distance = np.zeros(centerline_unsmoothed.shape[2]-1)
    distance_phys = np.zeros(centerline_unsmoothed.shape[2]-1)
    z_max = centerline_unsmoothed.shape[2]-1
    for i in range(z_max/2, z_max):

        ind1 = np.argmax(centerline_unsmoothed[:, :, i])
        X1,Y1 = ind2sub(centerline_unsmoothed[:, :, i].shape,ind1)
        coord_i_pix = [[X1, Y1, i]]
        coord_i = img_centerline_unsmoothed.transfo_pix2phys(coordi=coord_i_pix)

        ind2 = np.argmax(centerline_unsmoothed[:, :, i+1])
        X2,Y2 = ind2sub(centerline_unsmoothed[:, :, i+1].shape,ind2)
        coord_i_pix = [[X2, Y2, i+1]]
        coord_i_plus_1 = img_centerline_unsmoothed.transfo_pix2phys(coordi=coord_i_pix)

        distance_phys[i] = np.linalg.norm(np.array([coord_i_plus_1[0][0] - coord_i[0][0], coord_i_plus_1[0][1] - coord_i[0][1]]))
        distance[i] = np.linalg.norm(np.array([X2 - X1, Y2-Y1]))
    for i in range((z_max/2), 0, -1):

        ind1 = np.argmax(centerline_unsmoothed[:, :, i])
        X1,Y1 = ind2sub(centerline_unsmoothed[:, :, i].shape,ind1)
        coord_i_pix = [[X1, Y1, i]]
        coord_i = img_centerline_unsmoothed.transfo_pix2phys(coordi=coord_i_pix)

        ind2 = np.argmax(centerline_unsmoothed[:, :, i-1])
        X2,Y2 = ind2sub(centerline_unsmoothed[:, :, i-1].shape,ind2)
        coord_i_pix = [[X2, Y2, i-1]]
        coord_i_minus_1 = img_centerline_unsmoothed.transfo_pix2phys(coordi=coord_i_pix)

        distance_phys[i] = np.linalg.norm(np.array([coord_i_minus_1[0][0] - coord_i[0][0], coord_i_minus_1[0][1] - coord_i[0][1]]))
        distance[i] = np.linalg.norm(np.array([X2 - X1, Y2-Y1]))

    from scipy.stats import norm
    import matplotlib.pyplot as plt

    if std is None and mu is None:
        mu, std = norm.fit(distance)

    # removed_outlier = False

    for i in range(z_max/2, z_max):
        if distance_phys[i] > 10:
            ind1 = np.argmax(centerline_unsmoothed[:, :, i])
            X1,Y1 = ind2sub(centerline_unsmoothed[:, :, i].shape,ind1)
            ind2 = np.argmax(centerline_unsmoothed[:, :, i+1])
            X2,Y2 = ind2sub(centerline_unsmoothed[:, :, i+1].shape,ind2)
            centerline_unsmoothed[:,:,i+1] = 0
            # if np.amax(centerline_unsmoothed[:, :, i+1]) == 0:
            #     centerline_unsmoothed[X1,Y1,i+1] = 1
            # removed_outlier = True

    for i in range((z_max/2)-1, 0, -1):
        if distance_phys[i] > 10:
            ind1 = np.argmax(centerline_unsmoothed[:, :, i])
            X1,Y1 = ind2sub(centerline_unsmoothed[:, :, i].shape,ind1)
            ind2 = np.argmax(centerline_unsmoothed[:, :, i-1])
            X2,Y2 = ind2sub(centerline_unsmoothed[:, :, i-1].shape,ind2)
            centerline_unsmoothed[:, :, i-1] = 0
            # if np.amax(centerline_unsmoothed[:, :, i-1]) == 0:
            #     centerline_unsmoothed[X1, Y1, i-1] = 1
            # removed_outlier = True

    # if removed_outlier:
    #     img_centerline_unsmoothed.data = centerline_unsmoothed
    #     centerline_unsmoothed = remove_outliers(img_centerline_unsmoothed, std, mu)

    return centerline_unsmoothed

    # from scipy.stats import norm
    # import matplotlib.pyplot as plt
    # size = np.size(dx)
    # norms = np.zeros(size-1)
    # for i in range(0, np.size(dx)-1):
    #     array = np.array([dx[i], dy[i]])
    #     norms[i] = np.linalg.norm(array)
    #
    # mu, std = norm.fit(norms)
    # # plt.hist(norms, bins=0.1, normed=True, alpha=0.6, color='g')
    # return norms


def symmetry_detector_right_left(data, dim, y_crop_elevation=0, cropped_xy=0):
    """
    This function computes the cross correlation for L-R body symmetry for an RPI image
    This function expects an RPI image
    :param data: input image data
    :param dim: dimension of image
    :param y_crop_elevation: if the SC is not included when cropping on A-P, then set this parameter to include it. This parameter value is in pixel # TODO:Give a better explanation
    :return: L-R body symmetry (normalised between 0 and 1)
    """
    from scipy.ndimage.filters import gaussian_filter
    img_data = gaussian_filter(data, [0, 5, 5])

    if cropped_xy:
        x_mid = np.round(dim[0]/2)
        x_crop_min = int(x_mid - (0.35/2)*dim[0])
        x_crop_max = int(x_mid + (0.35/2)*dim[0])

        img_data[0:x_crop_min,:,:] = 0
        img_data[x_crop_max:-1,:,:] = 0

    slice_p = np.squeeze(np.sum(img_data, 1))
    # slice_p = np.squeeze(np.sum(img_data[x_crop_min:x_crop_max,:,:], 1)) # test
    slice_p_reversed = np.flipud(slice_p)
    m, n = slice_p.shape
    cross_corr = ((2*m)-1, n)
    cross_corr = np.zeros(cross_corr)
    for iz in range(0, np.size(slice_p[1])):
        corr1 = slice_p[:, iz]
        corr2 = slice_p_reversed[:, iz]
        cross_corr[:, iz] = np.double(np.correlate(corr1, corr2, "full"))
        max_value = np.max(cross_corr[:, iz])
        if max_value == 0:
            cross_corr[:, iz] = 0
        else:
            cross_corr[:, iz] = cross_corr[:, iz]/max_value
    data_out = np.zeros((dim[0], dim[2]))
    index1 = np.round(np.linspace(0,2*m-3, m))
    index2 = np.round(np.linspace(1,2*m-2, m))
    for i in range(0,m):
        indx1 = int(index1[i])
        indx2 = int(index2[i])
        out1 = cross_corr[indx1, :]
        out2 = cross_corr[indx2, :]
        data_out[i, :] = 0.5*(out1 + out2)
    result = np.hstack([data_out[:, np.newaxis, :] for i in range(0, dim[1])])

    return result


def symmetry_detector_anterior_posterior(data, dim):
    slice_p = np.squeeze(np.sum(data, 0))
    slice_p_reversed = np.flipud(slice_p)
    m, n= slice_p.shape
    cross_corr = ((2*m)-1, n)
    cross_corr = np.zeros(cross_corr)
    for iz in range(0, np.size(slice_p[1])):
        corr1 = slice_p[:, iz]
        corr2 = slice_p_reversed[:, iz]
        cross_corr[:, iz] = np.double(np.correlate(corr1, corr2, "full"))
        max_value = np.max(cross_corr[:, iz])
        cross_corr[:, iz] = cross_corr[:, iz]/max_value
    data_out = np.zeros((m, n))
    index1 = np.round(np.linspace(0,2*m-3, m))
    index2 = np.round(np.linspace(1,2*m-2, m))
    for i in range(0,m):
        indx1 = int(index1[i])
        indx2 = int(index2[i])
        out1 = cross_corr[indx1, :]
        out2 = cross_corr[indx2, :]
        data_out[i, :] = 0.5*(out1 + out2)
    result = np.hstack([data_out[:, np.newaxis, :] for i in range(0, dim[1])])
    return result


def ind2sub(array_shape, ind):
    rows = (ind.astype('int') / array_shape[1])
    cols = (ind.astype('int') % array_shape[1]) # or numpy.mod(ind.astype('int'), array_shape[1])
    return (rows, cols)


def get_centerline(data, dim):

    from scipy.ndimage.filters import gaussian_filter
    # if img is not Image:
    #     raise Exception("Input for get_centerline is not of type Image")
    centerline = np.zeros(dim)

    data = gaussian_filter(data, [1,1,0])

    for iz in range(0, dim[2]):
        ind = np.argmax(data[:, :, iz])
        X,Y = ind2sub(data[:, :, iz].shape,ind)
        centerline[X,Y,iz] = 1

        if iz == 1:
            centerline[X,Y,0] = 1

        if iz == dim[2]-1:
            centerline[:,:,dim[2]-1] = centerline[:, :, dim[2]-2]

    return centerline


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
    def __init__(self, input_image, contrast=None, verbose=1):
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


    def test_debug(self):
        import matplotlib.pyplot as plt
        path_tmp = sct.tmp_create()
        sct.tmp_copy_nifti(self.input_image.absolutepath, path_tmp, 'raw.nii')
        sct.run('cp imageVesselNessFilter.nii.gz '+path_tmp+'imageVesselNessFilter.nii.gz')
        os.chdir(path_tmp)

        img = Image("raw.nii")
        raw_orientation = img.change_orientation()
        block_sym = block_symmetry(img, img.dim[1]/3)
        img.data = symmetry_detector_right_left(img.data, img.dim)
        img.change_orientation(raw_orientation)
        img.file_name += "_symetry"
        img.save()
        raw_symetry_data = img.data
        img.data = block_sym
        img.change_orientation(raw_orientation)
        img.file_name = "block_symmetry"
        img.save()


        # img.data = ap_sym
        # img.change_orientation(raw_orientation)
        # img.file_name += "_ap_symetry"
        # img.save()

        img = Image('imageVesselNessFilter.nii.gz')
        raw_orientation = img.change_orientation()
        minimum_path_data = get_minimum_path(img.data, invert=1)
        img.data = minimum_path_data
        img.change_orientation()
        img.file_name = "min_path"
        img.save()
        # treshold_minimum_path = treshold_min_path(img.data, 0.95)
        normalised_symmetry = equalize_array_histogram(block_sym)
        # enhanced_symmetry = np.power(normalised_symmetry, 2)
        spine_detect_data = np.multiply(minimum_path_data, normalised_symmetry)
        img.data = spine_detect_data
        img.change_orientation('RPI')
        img.file_name = "symetry_weighted_min_path"
        img.save()

        centerline_unsmoothed = get_centerline(spine_detect_data, spine_detect_data.shape)
        img.data = centerline_unsmoothed
        img.change_orientation('RPI')
        img.file_name = "centerline_unsmoothed_outliers"
        img.save()
        img.data = remove_outliers(img)
        # norms = remove_outliers(img.data)
        img.change_orientation('RPI')
        img.file_name = "centerline_unsmoothed"
        img.save()

        x, y, z, dx, dy, dz = smooth_centerline("centerline_unsmoothed_outliers.nii.gz")
        # x, y, z, dx, dy, dz = smooth_centerline("centerline_unsmoothed.nii.gz")

        centerline_dim = img.dim

        img.data = np.zeros(centerline_dim)
        for i in range(0, np.size(x)-1):
            # if i == 0:
            #     img.data[int(x[i+1] - dx[i+1]), int(y[i+1]-dy[i+1]), int(z[i])] = 1
            # if i == np.size(x)-1:
            #     img.data[int(x[i-1] + dx[i-1]), int(y[i-1]+dy[i-1]), int(z[i])] = 1
            img.data[int(x[i]), int(y[i]), int(z[i])] = 1




        img.change_orientation(raw_orientation)
        img.file_name = "centerline"
        img.save()

        # copy back centerline
        os.chdir('../')
        # sct.tmp_copy_nifti(path_tmp + 'centerline.nii.gz',self.input_image.path,self.input_image.file_name+'_centerline'+self.input_image.ext)


    def execute(self):
        # TODO : Use symmetry detector class (integrate the function into the detector class)

        print 'Execution of the SCAD algorithm'

        if self.verbose:
            import matplotlib.pyplot as plt # import for debug purposes

        # create tmp and copy input
        path_tmp = sct.tmp_create()
        sct.tmp_copy_nifti(self.input_image.absolutepath, path_tmp, 'raw.nii')
        os.chdir(path_tmp)

        img = Image("raw.nii")
        raw_orientation = img.change_orientation()
        img.data = symmetry_detector_right_left(img.data, img.dim)
        img.change_orientation(raw_orientation)

        if self.verbose:
            img.file_name += "_symetry"
            img.save()

        # Save symmetry data in memory
        raw_symetry_data = img.data

        # vesselness filter
        sct.run('sct_vesselness -i raw.nii -t ' + self._contrast)

        # load vesselness filter data and perform minimum path on it
        img = Image('imageVesselNessFilter.nii.gz')
        raw_orientation = img.change_orientation()
        minimum_path_data = get_minimum_path(img.data, invert=1)

        if self.verbose:
            img.data = minimum_path_data
            img.change_orientation(raw_orientation)
            img.file_name = "min_path"
            img.save()

        # multiply normalised symmetry data with the minimum path result
        spine_detect_data = np.multiply(minimum_path_data, equalize_array_histogram(raw_symetry_data))
        if self.verbose:
            img.data = spine_detect_data
            img.change_orientation('RPI')
            img.file_name = "symetry_weighted_min_path"
            img.save()

        # symmetry_iter = symmetry_detector_right_left(spine_detect_data, spine_detect_data.shape) # test
        # if self.verbose:
        #     img.data = symmetry_iter
        #     img.change_orientation()
        #     img.file_name = "symmetry_iter"
        #     img.save()
        # get the centerline from the minimal path
        img.data = get_centerline(spine_detect_data, spine_detect_data.shape)
        img.file_name = "centerline_unsmoothed"
        img.save()

        # spine_detect_data = np.multiply(minimum_path_data, equalize_array_histogram(symmetry_iter))
        # img.data = get_centerline(spine_detect_data, spine_detect_data.shape)
        # img.file_name = "centerline_symmetry_iter"
        # img.save()

        x, y, z, dx, dy, dz = smooth_centerline("centerline_unsmoothed.nii.gz")

        centerline_dim = img.dim
        img.data = np.zeros(centerline_dim)
        for i in range(0, np.size(x)-1):
            img.data[int(x[i]), int(y[i]), int(z[i])] = 1

        img.change_orientation(raw_orientation)
        img.file_name = "centerline"
        img.save()

        # copy back centerline
        os.chdir('../')
        sct.tmp_copy_nifti(path_tmp + 'centerline.nii.gz',self.input_image.path,self.input_image.file_name+'_centerline'+self.input_image.ext)


class SymmetryDetector(Algorithm):
    def __init__(self, input_image, contrast=None, verbose=0, direction="lr"):
        super(SymmetryDetector, self).__init__(input_image)
        self._contrast = contrast
        self._verbose = verbose
        self.direction = direction

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
        """
        This function finds the line by line body symmetry from an input image
        :return: normalised cross-correlation values for the image
        """
        print 'Execution of the symmetry detection algorithm'



if __name__ == "__main__":
    parser = Script.get_parser()

    arguments = parser.parse(sys.argv[1:])

    input_image = Image(arguments["-i"])
    contrast_type = arguments["-t"]

    scad = SCAD(input_image)
    scad.contrast = contrast_type
    scad.execute()
