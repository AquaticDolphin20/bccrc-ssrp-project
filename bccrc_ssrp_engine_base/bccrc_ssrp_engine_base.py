# bccrc-ssrp engine base, by: Ian Janzen, 2025-06-10

import os
import pandas as pd
from pathlib import Path
import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt

import pydicom as pyd 


def add_element(dict_in: dict, key_in, value):
    if key_in not in dict_in:
        dict_in[key_in] = []
    dict_in[key_in].append(value)


def walkDir(path_in:str):
    # lists all filepaths
    list_out = []
    for dirpath, dirnames, filenames in os.walk(path_in):
        for filename in filenames:
                #  print(dirpath, filename)
                fp = os.path.join(dirpath, filename)
                list_out.append(fp)
    return list_out

def coord_conversion(val, position, pixel_sp):
    return int(distance.euclidean(val, position) / pixel_sp)


path_to_csv = r'..\\path\\to\\csv\\..'

df_cirrus = pd.read_csv(path_to_csv)



for ix in range(len(df_cirrus)): 
    
    # if you prefer a dict
    nodule_dict = {}
    _ = [add_element(nodule_dict, col, df_cirrus[col].iloc[ix]) for col in list(df_cirrus.columns)]
    pt = nodule_dict['PC_NUMBER']
    ndr = nodule_dict['NoduleID']
    px = float(nodule_dict['PixelSpacing'].split('\\')[-1])
    sl_thicc = float(nodule_dict['SliceThickness'])

    # if you prefer to pull the information you need
    pt = df_cirrus['PC_NUMBER'].iloc[ix] 
    ndr = df_cirrus['NoduleID'].iloc[ix]
    px = float(df_cirrus['PixelSpacing'].iloc[ix].split('\\')[-1])
    sl_thicc = float(df_cirrus['SliceThickness'].iloc[ix])

    print(f'Working on Patient#{pt} for nodule #{ndr}...')

    # kernel type:
    kern = df_cirrus['ConvolutionKernel'].iloc[ix]

    xc = df_cirrus.CoordX.iloc[ix]
    yc = df_cirrus.CoordY.iloc[ix]
    zc = df_cirrus.CoordZ.iloc[ix]

    # get the minor and major axis length
    min_axis = df_cirrus.MinorAxis_mm.iloc[ix] / 2 # radius
    max_axis = df_cirrus.MajorAxis_mm.iloc[ix] / 2 

    BBOX_Px_Z_CENTER = df_cirrus['SliceNumber'].iloc[ix] + 1 # ASSUMPTION: Cirrus uses 0-indexing, hence the +1
    BBOX_Z = [np.floor(BBOX_Px_Z_CENTER - (max_axis / px)), np.ceil(BBOX_Px_Z_CENTER + (max_axis / px))]

    # this method is not efficient as it will reload all the dicoms for each nodule, despite being from the same patient
    path_to_dicom = rf'..\\path\\to\\{pt}\\dicom\\'
    dicom_list = walkDir(path_to_dicom)

    # iterate through dicoms, find what convolution kernel you need
    volume = 0
    sl_last = 0
    for fp in dicom_list:
        # try to read dicom
        try:
            dcm = pyd.dcmread(fp)
        except pyd.InvalidDicomError as invDError:
        #   print(invDError)
            continue

        # try to read convolution kernel and image type
        try:
            conv_kernel = dcm.ConvolutionKernel  # this "filter" is a flawed process. See comment in error.
            img_type = str(dcm.ImageType)
        except AttributeError as e:
            print(e)  # This should not trigger by the time we're done as the BBOX might be need to be laid over CXR too

        # only keeps images with matching convolution kernel
        if conv_kernel != kern:
            continue

        dcm_im = dcm.pixel_array # should be the dicom image, likley in Hounsfield Units (HU)
        # it may also be the 'Pixel Data' tag
        # if you want to visualize the image, you will have to add the RescaleIntercept value to the image
        # and use the Lung window value, default is [-600, 1500] I believe, and then finally convert it to uint8 [px val 0-255]

        if not volume:
            im_shape = np.shape(dcm_im)
            volume = np.zeros((im_shape[0],im_shape[1],2000))
        
        dcm_keep = dcm
        PATIENT_POSITION = dcm['ImagePositionPatient'].value
        sl = dcm.InstanceNumber - 1 # or sometimes it is dcm.SliceNumber. Can't remember if this is 0-indexed or not, I'm assuming no, hence the -1.
        volume[:,:,sl] = dcm_im[:,:,:] # this might give you some headaches. You could vectorize it all then reshape back too.
        
        if sl > sl_last:
            sl_last = sl # keeping the last slice_value

        xc_conv = coord_conversion(xc, PATIENT_POSITION[0], px)
        yc_conv = coord_conversion(yc, PATIENT_POSITION[1], px)
        zc_conv = coord_conversion(zc, PATIENT_POSITION[2], sl_thicc)

        # I'm guessing here with this
        if BBOX_Z[0] <= zc_conv <= BBOX_Z[1]:
            do_something = 1
            #here we would maybe have to include the Graphics Annotation Tag for this particular slice

    volume = volume[:,:,0:sl_last] # there is a more elegant way of trimming the excess zeroes off the volume for viewing, but this gets the job done. 
    


