__author__ = 'taduv_admin'

from msct_image import Image
import numpy as np
import sct_straighten_spinalcord



def vertebral_detection(fname, fname_centerline, verbose=0):

    shift_AP = 14                                       # shift the centerline on the spine in mm default : 17 mm
    size_AP = 6                                         # mean around the centerline in the anterior-posterior direction in mm
    size_RL = 5                                         # mean around the centerline in the right-left direction in mm


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
    x, y, z, Tx, Ty, Tz = sct_straighten_spinalcord.smooth_centerline(fname_centerline)

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
    Wn = 0.02/frequency
    N = 2              #Order of the filter
    #    b, a = scipy.signal.butter(N, Wn, btype='low', analog=False, output='ba')
    b, a = scipy.signal.iirfilter(N, Wn, rp=None, rs=None, btype='high', analog=False, ftype='bessel', output='ba')
    I_detrend = scipy.signal.filtfilt(b, a, I[:,0], axis=-1, padtype='constant', padlen=None)
    I_detrend = I_detrend/(np.amax(I_detrend))
    I_detrend2= np.diff(I_detrend)



    #==================================================
    # step 1 : Find the First Peak
    #==================================================
    #locs=scipy.signal.find_peaks_cwt(np.squeeze(I_detrend2),np.arange(1,I_detrend2.size),gap_thresh=0.1)
    I_detrend2[I_detrend2<0.2]=0
    locs=np.squeeze(scipy.signal.argrelextrema(I_detrend2,np.greater))


    # remove peaks that are too closed
    locsdiff=np.diff(locs)
    ind=locsdiff>10
    locs = np.hstack((locs[ind],locs[-1]))

    if verbose:
        plt.figure()
        plt.plot(I_detrend2)
        plt.plot(locs,I_detrend2[locs],'+')

    # mean_distance_dict = scipy.io.loadmat('/home/django/kraju/code/spinalcordtoolbox_dev/src/vertebral_labeling/mean_distance.mat')
    # mean_distance = (mean_distance_dict.values()[2]).T
    # C1C2_distance = mean_distance[0:2]
    # mean_distance = mean_distance[level_start-1:len(mean_distance)-1]
    #
    # space = np.linspace(-5/scales[2], 5/scales[2], round(11/scales[2]), endpoint=True)
    # pattern = (np.sinc((space*scales[2])/15))**(20)
    # xmax_pattern = np.argmax(pattern)


    #correlation between the pattern and intensity profile
    #corr_all = scipy.signal.correlate(pattern,I_detrend[:,0])
    #corr_all = matplotlib.pyplot.xcorr(pattern,I_detrend[:,0])

    # pattern1 =  np.concatenate((pattern,np.zeros(len(I_detrend[:,0])-len(pattern))))
    # corr_all = scipy.signal.correlate(I_detrend[:,0],pattern1)
    # loc_corr = np.arange(-np.round((len(corr_all)/2)),np.round(len(corr_all)/2)+2)
    # index_fp = 0
    # count = 0
    # for i in range(len(corr_all)):
    #     if corr_all[i]>0.1:
    #         if i==0:
    #             if corr_all[i]<corr_all[i+1]:
    #                 index_fp = i
    #                 count = count + 1
    #         elif i==(len(corr_all)-1):
    #             if corr_all[i]<corr_all[i-1]:
    #                 index_fp = np.resize(index_fp,count+1)
    #                 index_fp[len(index_fp)-1] = i
    #         else:
    #             if corr_all[i]<corr_all[i+1]:
    #                 index_fp = np.resize(index_fp,count+1)
    #                 index_fp[len(index_fp)-1] = i
    #                 count = count + 1
    #             elif corr_all[i]<corr_all[i-1]:
    #                 index_fp = np.resize(index_fp,count+1)
    #                 index_fp[len(index_fp)-1] = i
    #                 count = count + 1
    #     else:
    #         if i==0:
    #             index_fp = i
    #             count = count + 1
    #         else:
    #             index_fp = np.resize(index_fp,count+1)
    #             index_fp[len(index_fp)-1] = i
    #             count = count + 1
    #
    #
    # mask_fp = np.ones(len(corr_all), dtype=bool)
    # mask_fp[index_fp] = False
    # value = corr_all[mask_fp]
    # loc_corr = loc_corr[mask_fp]
    #
    # loc_corr = loc_corr - I_detrend.shape[0]
    # loc_first_peak = xmax_pattern - loc_corr[np.amax(np.where(value>1))]
    # Mcorr1 = value[np.amax(np.where(value>1))]
    #
    # #building the pattern that has to be added at each iteration in step 2
    #
    # if xmax_pattern<loc_first_peak:
    #     template_truncated = np.concatenate((np.zeros((loc_first_peak-xmax_pattern)),pattern))
    #
    # else:
    #     template_truncated = pattern[(xmax_pattern-loc_first_peak-1):]
    # xend = np.amax(np.where(template_truncated>0.02))
    # pixend = xend - loc_first_peak
    #
    # if label.verbose==1:
    #     pl.plot(template_truncated)
    #     pl.plot(I_detrend)
    #     pl.title('Detection of First Peak')
    #     pl.xlabel('direction anterior-posterior (mm)')
    #     pl.ylabel('intensity')
    #     pl.show()
    #
    # loc_peak_I = np.arange(len(I_detrend[:,0]))
    # count = 0
    # index_p = 0
    # for i in range(len(I_detrend[:,0])):
    #     if I_detrend[i]>0.15:
    #         if i==0:
    #             if I_detrend[i,0]<I_detrend[i+1,0]:
    #                 index_p = i
    #                 count  =  count + 1
    #         elif i==(len(I_detrend[:,0])-1):
    #             if I_detrend[i,0]<I_detrend[i-1,0]:
    #                 index_p = np.resize(index_p,count+1)
    #                 index_p[len(index_p)-1] = i
    #         else:
    #             if I_detrend[i,0]<I_detrend[i+1,0]:
    #                 index_p = np.resize(index_p,count+1)
    #                 index_p[len(index_p)-1] = i
    #                 count = count+1
    #             elif I_detrend[i,0]<I_detrend[i-1,0]:
    #                 index_p = np.resize(index_p,count+1)
    #                 index_p[len(index_p)-1] = i
    #                 count = count+1
    #     else:
    #         if i==0:
    #             index_p = i
    #             count  =  count + 1
    #         else:
    #             index_p = np.resize(index_p,count+1)
    #             index_p[len(index_p)-1] = i
    #             count = count+1
    #
    # mask_p = np.ones(len(I_detrend[:,0]), dtype=bool)
    # mask_p[index_p] = False
    # value_I = I_detrend[mask_p]
    # loc_peak_I = loc_peak_I[mask_p]
    #
    # count = 0
    # for i in range(len(loc_peak_I)-1):
    #     if i==0:
    #         if loc_peak_I[i+1]-loc_peak_I[i]<round(10/scales[1]):
    #             index = i
    #             count = count + 1
    #     else:
    #         if (loc_peak_I[i+1]-loc_peak_I[i])<round(10/scales[1]):
    #             index =  np.resize(index,count+1)
    #             index[len(index)-1] = i
    #             count = count + 1
    #         elif (loc_peak_I[i]-loc_peak_I[i-1])<round(10/scales[1]):
    #             index =  np.resize(index,count+1)
    #             index[len(index)-1] = i
    #             count = count + 1
    #
    # mask_I = np.ones(len(value_I), dtype=bool)
    # mask_I[index] = False
    # value_I = value_I[mask_I]
    # loc_peak_I = loc_peak_I[mask_I]
    #
    # from scipy.interpolate import UnivariateSpline
    # fit = UnivariateSpline(loc_peak_I,value_I)
    # P = fit(np.arange(len(I_detrend)))
    #
    # for i in range(len(I_detrend)):
    #     if P[i]>0.1:
    #         I_detrend[i,0] = I_detrend[i,0]/P[i]
    #
    # if label.verbose==1:
    #     pl.xlim(0,len(I_detrend)-1)
    #     pl.plot(loc_peak_I,value_I)
    #     pl.plot(I_detrend)
    #     pl.plot(P,color='y')
    #     pl.title('Setting values of peaks at one by fitting a smoothing spline')
    #     pl.xlabel('direction superior-inferior (mm)')
    #     pl.ylabel('normalized intensity')
    #     pl.show(block=False)
    #
    # #===================================================================================
    # # step 2 : Cross correlation between the adjusted template and the intensity profile
    # #          local moving of template's peak from the first peak already found
    # #===================================================================================
    #
    # mean_distance_new = mean_distance
    # mean_ratio = np.zeros(len(mean_distance))
    # L = np.round(1.2*max(mean_distance)) - np.round(0.8*min(mean_distance))
    # corr_peak  = np.zeros((L,len(mean_distance)))
    #
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
    #         corr_template = scipy.signal.correlate(I_detrend[:,0],template_resize_peak)
    #
    #         if len(I_detrend[:,0])>len(template_resize_peak):
    #             val = np.dot(I_detrend[:,0],template_resize_peak1.T)
    #         else:
    #             I_detrend_2 = np.concatenate((I_detrend[:,0],np.zeros(len(template_resize_peak)-len(I_detrend[:,0]))))
    #             val = np.dot(I_detrend_2,template_resize_peak.T)
    #         corr_peak[i_scale,i_peak] = val
    #
    #         if label.verbose==1:
    #             pl.xlim(0,len(I_detrend[:,0]))
    #             pl.plot(I_detrend[:,0])
    #             pl.plot(template_resize_peak)
    #             pl.show(block=False)
    #
    #             pl.plot(corr_peak[:,i_peak],marker='+',linestyle='None',color='r')
    #             pl.title('correlation value against the displacement of the peak (px)')
    #             pl.show(block=False)
    #
    #     max_peak = np.amax(corr_peak[:,i_peak])
    #     index_scale_peak = np.where(corr_peak[:,i_peak]==max_peak)
    #     good_scale_peak = scale_peak[index_scale_peak][0]
    #     Mcorr = Mcorr1
    #     Mcorr = np.resize(Mcorr,i_peak+2)
    #     Mcorr[i_peak+1] = np.amax(corr_peak[:,0:(i_peak+1)])
    #     flag = 0
    #
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
    #     mean_distance_new[i_peak] = good_scale_peak + xmax_pattern + pixend
    #     mean_ratio[i_peak] = np.mean(mean_distance_new[:,0:i_peak]/mean_distance[:,0:i_peak])
    #
    #     template_truncated = template_resize_peak
    #
    #     if label.verbose==1:
    #         pl.plot(I_detrend[:,0])
    #         pl.plot(template_truncated)
    #         pl.xlim(0,(len(I_detrend[:,0])-1))
    #         pl.show()
    #
    # minpeakvalue = 0.5
    # loc_disk = np.arange(len(template_truncated))
    # count = 0
    # index_disk = 0
    # for i in range(len(template_truncated)):
    #     if template_truncated[i]>=minpeakvalue:
    #         if i==0:
    #             if template_truncated[i]<template_truncated[i+1]:
    #                 index_disk = i
    #                 count  =  count + 1
    #         elif i==(len(template_truncated)-1):
    #             if template_truncated[i]<template_truncated[i-1]:
    #                 index_disk = np.resize(index_disk,count+1)
    #                 index_disk[len(index_disk)-1] = i
    #         else:
    #             if template_truncated[i]<template_truncated[i+1]:
    #                 index_disk = np.resize(index_disk,count+1)
    #                 index_disk[len(index_disk)-1] = i
    #                 count = count+1
    #             elif template_truncated[i]<template_truncated[i-1]:
    #                 index_disk = np.resize(index_disk,count+1)
    #                 index_disk[len(index_disk)-1] = i
    #                 count = count+1
    #     else:
    #         if i==0:
    #             index_disk = i
    #             count  =  count + 1
    #         else:
    #             index_disk = np.resize(index_disk,count+1)
    #             index_disk[len(index_disk)-1] = i
    #             count = count+1
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
        for iz in range(locs[i],locs[i+1]):
            centerline.data[np.round(x[iz]),np.round(y[iz]),iz]=i+2
    for iz in range(locs[-1],len(z)):
            centerline.data[np.round(x[iz]),np.round(y[iz]),iz]=i+3

    #centerline.change_orientation(raw_orientation)
    centerline.file_name+='_labeled'
    centerline.save()

    return locs