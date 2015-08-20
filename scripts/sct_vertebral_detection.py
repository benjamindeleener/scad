#!/usr/bin/env python
# __author__ = 'taduv_admin'

from msct_base_classes import BaseScript
from msct_image import Image
import numpy as np
import sct_straighten_spinalcord
import sct_utils as sct
from msct_parser import Parser
import sys


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
        parser.add_option(name="-centerline",
                          type_value="file",
                          description="Binary image of the centerline",
                          mandatory=True,
                          example="centerline.nii.gz")
        parser.add_option(name="-seg",
                          type_value="file",
                          description="input image.",
                          mandatory=True,
                          example="segmentation.nii.gz")
        parser.add_option(name="-v",
                          type_value="multiple_choice",
                          description="1: display on, 0: display off (default)",
                          mandatory=False,
                          example=["0", "1"],
                          default_value="1")
        parser.add_option(name="-t",
                          type_value="multiple_choice",
                          description="Image Contrast",
                          mandatory=True,
                          example=["T1", "T2"])
        parser.add_option(name="-h",
                          type_value=None,
                          description="display this help",
                          mandatory=False)
        return parser


def vertebral_detection(fname, fname_centerline, fname_segmentation=None, contrast='T1', verbose=0):

    shift_AP = 14                                       # shift the centerline on the spine in mm default : 17 mm
    size_AP = 3                                         # mean around the centerline in the anterior-posterior direction in mm
    size_RL = 3                                         # mean around the centerline in the right-left direction in mm


    if verbose:
        import matplotlib.pyplot as plt

    img = Image(fname)
    img.change_orientation()
    # # extract path/file/extension
    # path_anat, file_anat, ext_anat = sct.extract_fname(input_anat)
    # path_centerline, file_centerline, ext_centerline = sct.extract_fname(input_centerline)


    #==================================================
    # Calculation of the profile intensity
    #==================================================

    shift_AP = shift_AP*img.pixdim[1]
    size_AP = size_AP*img.pixdim[1]
    size_RL = size_RL*img.pixdim[0]

    # normal vector of the orthogonal plane to the centerline i.e tangent vector to the centerline
    path_centerline, file_centerline, ext_centerline = sct.extract_fname(fname_centerline)
    sct.run('sct_orientation -i ' + fname_centerline + ' -s RPI')
    x, y, z, Tx, Ty, Tz = sct_straighten_spinalcord.smooth_centerline(path_centerline + file_centerline + '_RPI' + ext_centerline)

    #build intensity profile along the centerline
    I = np.zeros((len(y),1))

    #  mask where intensity profile will be taken
    if verbose:
        mat=img.copy()
        mat.data=np.zeros(mat.dim)


    for iz in range(len(z)):


        # vector perpendicular to the cord in RL direction
        P1 = np.array([1, 0, -Tx[iz]/Tz[iz]])
        P1 = P1/np.linalg.norm(P1)
        # vector perpendicular to the cord in AP direction
        P2 = np.array([0, 1, -Ty[iz]/Tz[iz]])
        P2 = P2/np.linalg.norm(P2)

        indexRL = range(-np.int(round(size_RL)),np.int(round(size_RL)))
        indexAP = range(0,np.int(round(size_AP)))+np.array(shift_AP)

        # loop over coordinates of perpendicular plane
        for i_RL in indexRL:
            for i_AP in indexAP:
                i_vect = np.round(np.array([x[iz],y[iz],z[iz]])+P1*i_RL+P2*i_AP)
                i_vect = np.minimum(np.maximum(i_vect,0),np.array(img.dim)-1) # check if index stays in image dimension
                I[iz] = I[iz] + img.data[i_vect[0],i_vect[1],i_vect[2]]

                # create a mask with this perpendicular plane
                if verbose:
                    mat.data[i_vect[0],i_vect[1],i_vect[2]] = 1

    if verbose:
        mat.file_name='mask'
        mat.save()

    # Detrending Intensity
    start_centerline_y = y[0]
    X = np.where(I==0)
    mask2 = np.ones((len(y),1), dtype=bool)
    mask2[X,0] = False

    # low pass filtering
    import scipy.signal
    frequency = 2/img.pixdim[2]
    Wn = 0.1/frequency
    N = 2              #Order of the filter
    #    b, a = scipy.signal.butter(N, Wn, btype='low', analog=False, output='ba')
    b, a = scipy.signal.iirfilter(N, Wn, rp=None, rs=None, btype='high', analog=False, ftype='bessel', output='ba')
    I_detrend = scipy.signal.filtfilt(b, a, I[:,0], axis=-1, padtype='constant', padlen=None)
    I_detrend = I_detrend/(np.amax(I_detrend))


    #==================================================
    # step 1 : Find the First Peak
    #==================================================
    if contrast == 'T1':
        I_detrend2= np.diff(I_detrend)
    elif contrast == 'T2':
        space = np.linspace(-10/img.pixdim[2], 10/img.pixdim[2], round(21/img.pixdim[2]), endpoint=True)
        pattern = (np.sinc((space*img.pixdim[2])/20))**(20)
        I_corr = scipy.signal.correlate(-I_detrend.squeeze().squeeze()+1,pattern,'same')
        b, a = scipy.signal.iirfilter(N, Wn, rp=None, rs=None, btype='high', analog=False, ftype='bessel', output='ba')
        I_detrend2 = scipy.signal.filtfilt(b, a, I_corr, axis=-1, padtype='constant', padlen=None)

    I_detrend2[I_detrend2<0.2]=0
    ind_locs=np.squeeze(scipy.signal.argrelextrema(I_detrend2,np.greater))



    # remove peaks that are too closed
    locsdiff=np.diff(z[ind_locs])
    ind=locsdiff>10
    ind_locs = np.hstack((ind_locs[ind],ind_locs[-1]))
    locs=z[ind_locs]

    if verbose:
        plt.figure()
        plt.plot(I_detrend2)
        plt.plot(ind_locs,I_detrend2[ind_locs],'+')



    #=====================================================================================
    # step 2 : Cross correlation between the adjusted template and the intensity profile.
    #          Local moving of template's peak from the first peak already found
    #=====================================================================================

    #For each loop, a peak is located at the most likely postion and then local adjustement is done.
    #The position of the next peak is calculated from previous positions

    sct.printv('\nFinding Cross correlation between the adjusted template and the intensity profile...',verbose)

    mean_distance = [12.1600, 20.8300, 18.0000, 16.0000, 15.1667, 15.3333, 15.8333,   18.1667,   18.6667,   18.6667,
    19.8333,   20.6667,   21.6667,   22.3333,   23.8333,   24.1667,   26.0000,   28.6667,   30.5000,   33.5000,
    33.0000,   31.3330]
    #
    # #Creating pattern
    space = np.linspace(-10/img.pixdim[2], 10/img.pixdim[2], round(21/img.pixdim[2]), endpoint=True)
    pattern = (np.sinc((space*img.pixdim[2])/20))**(20)
    I_corr = scipy.signal.correlate(I_detrend2.squeeze().squeeze()+1,pattern,'same')
    #
    # level_start=1
    # if contrast == 'T1':
    #     mean_distance = mean_distance[level_start-1:len(mean_distance)]
    #     xmax_pattern = np.argmax(pattern)
    # else:
    #     mean_distance = mean_distance[level_start+1:len(mean_distance)]
    #     xmax_pattern = np.argmin(pattern)          # position of the peak in the pattern
    # pixend = len(pattern) - xmax_pattern       #number of pixel after the peaks in the pattern
    #
    #
    # mean_distance_new = mean_distance
    # mean_ratio = np.zeros(len(mean_distance))
    #
    # L = np.round(1.2*max(mean_distance)) - np.round(0.8*min(mean_distance))
    # corr_peak  = np.zeros((L,len(mean_distance)))          # corr_peak  = np.nan #for T2
    #
    # #loop on each peak
    # for i_peak in range(len(mean_distance)):
    #     scale_min = np.round(0.80*mean_distance_new[i_peak]) - xmax_pattern - pixend
    #     if scale_min<0:
    #         scale_min = 0
    #
    #     scale_max = np.round(1.2*mean_distance_new[i_peak]) - xmax_pattern - pixend
    #     scale_peak = np.arange(scale_min,scale_max+1)
    #
    #     for i_scale in range(len(scale_peak)):
    #         template_resize_peak = np.concatenate([template_truncated,np.zeros(scale_peak[i_scale]),pattern])
    #         if len(I_detrend[:,0])>len(template_resize_peak):
    #             template_resize_peak1 = np.concatenate((template_resize_peak,np.zeros(len(I_detrend[:,0])-len(template_resize_peak))))
    #
    #         #cross correlation
    #         corr_template = scipy.signal.correlate(I_detrend[:,0],template_resize_peak)
    #
    #         if len(I_detrend[:,0])>len(template_resize_peak):
    #             val = np.dot(I_detrend[:,0],template_resize_peak1.T)
    #         else:
    #             I_detrend_2 = np.concatenate((I_detrend[:,0],np.zeros(len(template_resize_peak)-len(I_detrend[:,0]))))
    #             val = np.dot(I_detrend_2,template_resize_peak.T)
    #         corr_peak[i_scale,i_peak] = val
    #
    #         if verbose:
    #             plt.xlim(0,len(I_detrend[:,0]))
    #             plt.plot(I_detrend[:,0])
    #             plt.plot(template_resize_peak)
    #             plt.show(block=False)
    #
    #             plt.plot(corr_peak[:,i_peak],marker='+',linestyle='None',color='r')
    #             plt.title('correlation value against the displacement of the peak (px)')
    #             plt.show(block=False)
    #
    #     max_peak = np.amax(corr_peak[:,i_peak])
    #     index_scale_peak = np.where(corr_peak[:,i_peak]==max_peak)
    #     good_scale_peak = scale_peak[index_scale_peak][0]
    #     Mcorr = Mcorr1
    #     Mcorr = np.resize(Mcorr,i_peak+2)
    #     Mcorr[i_peak+1] = np.amax(corr_peak[:,0:(i_peak+1)])
    #     flag = 0
    #
    #     #If the correlation coefficient is too low, put the peak at the mean position
    #     if i_peak>0:
    #         if (Mcorr[i_peak+1]-Mcorr[i_peak])<0.4*np.mean(Mcorr[1:i_peak+2]-Mcorr[0:i_peak+1]):
    #             test = i_peak
    #             template_resize_peak = np.concatenate((template_truncated,np.zeros(round(mean_distance[i_peak])-xmax_pattern-pixend),pattern))
    #             good_scale_peak = np.round(mean_distance[i_peak]) - xmax_pattern - pixend
    #             flag = 1
    #     if i_peak==0:
    #         if (Mcorr[i_peak+1] - Mcorr[i_peak])<0.4*Mcorr[0]:
    #             template_resize_peak = np.concatenate((template_truncated,np.zeros(round(mean_distance[i_peak])-xmax_pattern-pixend),pattern))
    #             good_scale_peak = round(mean_distance[i_peak]) - xmax_pattern - pixend
    #             flag = 1
    #     if flag==0:
    #         template_resize_peak=np.concatenate((template_truncated,np.zeros(good_scale_peak),pattern))
    #
    #     #update mean-distance by a adjustement ratio
    #     mean_distance_new[i_peak] = good_scale_peak + xmax_pattern + pixend
    #     mean_ratio[i_peak] = np.mean(mean_distance_new[:,0:i_peak]/mean_distance[:,0:i_peak])
    #
    #     template_truncated = template_resize_peak
    #
    #     if verbose:
    #         plt.plot(I_detrend[:,0])
    #         plt.plot(template_truncated)
    #         plt.xlim(0,(len(I_detrend[:,0])-1))
    #         plt.show()
    #
    # #finding the maxima of the adjusted template
    # minpeakvalue = 0.5
    # loc_disk = np.arange(len(template_truncated))
    # index_disk = []
    # for i in range(len(template_truncated)):
    #     if template_truncated[i]>=minpeakvalue:
    #         if i==0:
    #             if template_truncated[i]<template_truncated[i+1]:
    #                 index_disk.append(i)
    #         elif i==(len(template_truncated)-1):
    #             if template_truncated[i]<template_truncated[i-1]:
    #                 index_disk.append(i)
    #         else:
    #             if template_truncated[i]<template_truncated[i+1]:
    #                 index_disk.append(i)
    #             elif template_truncated[i]<template_truncated[i-1]:
    #                 index_disk.append(i)
    #     else:
    #         index_disk.append(i)
    #
    # mask_disk = np.ones(len(template_truncated), dtype=bool)
    # mask_disk[index_disk] = False
    # loc_disk = loc_disk[mask_disk]
    # X1 = np.where(loc_disk > I_detrend.shape[0])
    # mask_disk1 = np.ones(len(loc_disk), dtype=bool)
    # mask_disk1[X1] = False
    # loc_disk = loc_disk[mask_disk1]
    # loc_disk = loc_disk + start_centerline_y - 1


    #=====================================================================
    # Step 3: Building of the labeled centerline and surface
    #=====================================================================

    # Project vertebral levels back to the centerline
    centerline = Image(fname_centerline)
    raw_orientation = centerline.change_orientation()
    centerline.data[:,:,:] = 0
    for iz in range(locs[0]):
            centerline.data[np.round(x[iz]),np.round(y[iz]),iz]=1
    for i in range(len(locs)-1):
        for iz in range(locs[i],min(locs[i+1],len(z))):
            centerline.data[np.round(x[iz]),np.round(y[iz]),iz]=i+2
    for iz in range(locs[-1],len(z)):
            centerline.data[np.round(x[iz]),np.round(y[iz]),iz]=i+3

    #centerline.change_orientation(raw_orientation)
    centerline.file_name+='_labeled'
    centerline.save()


    # color the segmentation with vertebral number
    if fname_segmentation:
        seg=Image(fname_segmentation)
        seg_raw_orientation = seg.change_orientation()
        x_seg,y_seg,z_seg=np.where(seg.data)
        for ivox in range(len(x_seg)):
            vox_coord = np.array([x_seg[ivox], y_seg[ivox], z_seg[ivox]])
            for iplane in range(len(locs)):
                ind=np.where(z==locs[iplane])
                vox_vector = vox_coord - np.hstack((x[ind],y[ind],z[ind]))
                normal2plane_vector = np.hstack((Tx[ind], Ty[ind], Tz[ind]))

                # if voxel is above the plane --> give the number of the plane
                if np.dot(vox_vector, normal2plane_vector) > 0:
                    seg.data[vox_coord[0], vox_coord[1], vox_coord[2]] = iplane+2

                else: # if the voxel gets below the plane --> next voxel
                    break
        seg.change_orientation(seg_raw_orientation)
        seg.file_name+='_labeled'
        seg.save()

    return locs

if __name__ == '__main__':
    parser = Script.get_parser()
    arguments = parser.parse(sys.argv[1:])

    vertebral_detection(arguments["-i"],arguments["-centerline"],fname_segmentation=arguments["-seg"],contrast=arguments["-t"], verbose=1)
