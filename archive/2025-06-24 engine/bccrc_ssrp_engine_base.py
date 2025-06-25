# bccrc-ssrp engine base, by: Ian Janzen, 2025-06-10

import os
import pandas as pd
from pathlib import Path
import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt

import pydicom as pyd 
from pydicom.misc import is_dicom
import highdicom as hdicom

# reads in all dicoms
def load_dicom(folder_path):
    dicom_datasets = [] 
    for filename in os.listdir(folder_path):
        filepath = os.path.join(folder_path, filename)

        # filter for DICOM files add on: 
        # if all dicom files end with .dcm, this is a useful filter to decrease runtime
        if not filename.lower().endswith('.dcm'): 
            continue

        # filter for DICOM file option 1
        try:
            ds = pydicom.dcmread(filepath) # can add argument stop_before_pixels = True to dcmread if we don't need to load all images
            dicom_datasets.append(ds)
        except pyd.InvalidDicomError as e: # if not dicom file, print error
            print("Invalid DICOM: ", e)
            continue

        # filter for DICOM file option 2
        if is_dicom(path): 
            ds = pydicom.dcmread(path)

    return dicom_datasets # returns datasets, not filepaths

dataset = load_dicom('test_data/52100009_baseline')

# reads in cirrus csv
path_to_csv = r'..\\path\\to\\csv\\..'

df_cirrus = pd.read_csv(path_to_csv)


# functions below are involved with bounding box creator and filtering process
# returns dicom images that has desired (can change input) dicom_attr match desired (can change input) df_attr from df_cirrus
def filter_images(dicom_folder, dicom_attr: str, desired_value: str, need_image: bool)
    # can just have df_attr be a straight up attribute of df


    # search dicom directory for patient dicom folder using PC_number
    filtered_dicom_datasets = [] # list will hold dicom folders matching desired dicom_attr
    for filename in os.listdir(folder_path):
        filepath = os.path.join(folder_path, filename)

        # filter for DICOM files add on: 
        # if all dicom files end with .dcm, this is a useful filter to decrease runtime
        if not filename.lower().endswith('.dcm'): 
            continue

        # filter for DICOM file option 1
        try:
            if need_image == True:
                ds = pydicom.dcmread(filepath)
            else: 
                ds = pydicom.dcmread(filepath, stop_before_pixels = True) # reduce run time by not loading images where image loading isn't necessary

            # filter process
            if getattr(ds, dicom_attr, None) == desired_value: # need to unstring dicom_attr for df_cirrus.dicom_attr
                patient_dicom_datasets.append(filepath)

        except pyd.InvalidDicomError as e: # if not dicom file, print error
            print("Invalid DICOM: ", e)
            continue

    return filtered_dicom_datasets

def find_slices(dicom_folder, df_cirrus, row)
    # takes in df_cirrus, required columns: slice, slice thickness, diameter..?, CoordZ..?
    # needs to convert from slice and slice thickness, which is pixel spacing, to patient position
    
    # finds which slices of dicom folder (organized image by axial slice) are necessary
    
    desired_value = df_cirrus.PC_NUMBER # we can change attribute later, but in this case we're using PC_Number of df_cirrus as identifier

    patient_dicom_datasets = filter_images(dicom_folder, "PatientID", desired_value, True) # whether PC_NUMBER is a str is prone to change

    # finds slices
    slices = []
    min_slice = np.floor(df.cirrus.SliceNumber - df.cirrus.Diameter_mm / df.SliceThickness)
    max_slice = np.ceil(df.cirrus.SliceNumber + df.cirrus.Diameter_mm / df.SliceThickness)
    for slice_num in range(min_slice, max_slice + 1):
    
    # this process only works if we assume all dicom folders will be named nicely
    dicom_slices = filter_images(patient_dicom_datasets, "")
        


    2. locate the slice number for patient dicom folder using slice number + instance
    
    df_cirrus.SliceNumber

    # finds which 

def create_bbox(df_cirrus, row) 
    # takes in df_cirrus, required columns: major axis, CoordX, CoordY, CoordZ, diameter
    # outputs a dictionary 
    
    bbox_dict = {} # for now will be dictionary of 3d coord, will make into numpy later

    min_axis = df_cirrus.MinorAxis_mm.iloc[row] / 2 # gets radius across minor axis - do we need this
    major_axis = df_cirrus.MajorAxis_mm.iloc[row] / 2
    bbox_z_center = df_cirrus['SliceNumber'].iloc[row] + 1 # +1 is there since cirrus uses 0-indexing.. do we even need this row
    
    # outputs a dictionary with keys x_coord, y_coord, z_coord each with a value of a set of [min, max] according to patient position
    bbox_dict[x_coord] = [np.floor(df_cirrus.CoordX.iloc[row] - (max_axis)), # values should be in patient position already, no need to divide by pixel spacing
                          np.ceil(df_cirrus.CoordX.iloc[row] + (max_axis))]
    bbox_dict[y_coord] = [np.floor(df_cirrus.CoordY.iloc[row] - (max_axis)),
                          np.ceil(df_cirrus.CoordY.iloc[row] + (max_axis))]
    bbox_dict[z_coord] = [df_cirrus.CoordZ.iloc[row] - df_cirrus.Diameter_mm.iloc[row], # do these rows need np.floor or np.ceil to round or is cirrus diameter accurate
                          df_cirrus.CoordZ.iloc[row] + df_cirrus.Diameter_mm.iloc[row]]

    return bbox_dict # is return statement necessary? probably not


# cirrus loop
for ix in range(len(df_cirrus)):



# match cirrus pancan with patient ID and create bounding box
    nodule_dict = {} # creates/resets dictionary
    add_element(nodule_dict, col, df_cirrus[col].iloc[ix]) for col in list(df_cirrus.columns) # adds all columns of dataset as keys and corresponding rows as values into dictionary
        # do we actually need to add to dictionary
    pc_num = nodule_dict['PC_NUMBER']
    ndr = nodule_dict['NoduleID']
    px = float(nodule_dict['PixelSpacing'].split('\\')[-1])
    sl_thicc = float(nodule_dict['SliceThickness'])

    create_bbox(df_cirrus, ix)
    matching_dicoms = search_dicom_folder(dataset, pc_num) # code function later
    image_filepath = locate_dicom_image(dataset, pc_num, slice_num) # code function later
    dicom_located_image = pydicom.dcmread(image_filepath) # code function later.. can I combine with previous function

# build_sr_file - does it require dicom image file (to reference metadata) or does it require UID?






highdicom.sr.VolumetricROIMeasurementsAndQualitativeEvaluations module needs to be used



graphic_data = np.array([
    (165.0, 200.0, 134.0),
    (170.0, 200.0, 134.0),
    (170.0, 220.0, 134.0),
    (165.0, 220.0, 134.0),
    (165.0, 200.0, 134.0),
])

# patient_id = "xxxxxx"
# image_id is encoded above
matches = [ds for ds in dicom_datasets if getattr(ds, "PatientID", None) == patient_id]

# data engine process: 
# 1. dicom loader (do we want to load one patient at a time or multiple patients)
# 2. csv loader for cirrus
# 3. coordinate conversion from cirrus pixel coordinates to patient orientation
# 4. checker that matches each csv row to respective dicom in patient
# 5. bounding box creator (maybe have a loop that checks over cirrus file to match patient dicom - so for one patient in dicom, do all bounding boxes of that patient that you find in cirrus first)
# 6. sr file creator to put dicom back in database


# def add_element(dict_in: dict, key_in, value):
#     if key_in not in dict_in:
#         dict_in[key_in] = []
#     dict_in[key_in].append(value)


# def walkDir(path_in:str):
#     # lists all filepaths
#     list_out = []
#     for dirpath, dirnames, filenames in os.walk(path_in):
#         for filename in filenames:
#                 #  print(dirpath, filename)
#                 fp = os.path.join(dirpath, filename)
#                 list_out.append(fp)
#     return list_out

# def coord_conversion(val, position, pixel_sp):
#     # calculates real world distance between val and pos
#     return int(distance.euclidean(val, position) / pixel_sp)


# path_to_csv = r'..\\path\\to\\csv\\..'

# df_cirrus = pd.read_csv(path_to_csv)



# for ix in range(len(df_cirrus)): # iterates through each row in df_cirrus
    
#     # if you prefer a dict
#     nodule_dict = {} # creates/resets dictionary
#     add_element(nodule_dict, col, df_cirrus[col].iloc[ix]) for col in list(df_cirrus.columns) # adds all columns of dataset as keys and corresponding value of row into dictionary
#     pt = nodule_dict['PC_NUMBER'] # assigns patient number of that row to pt
#     ndr = nodule_dict['NoduleID']
#     px = float(nodule_dict['PixelSpacing'].split('\\')[-1])
#     sl_thicc = float(nodule_dict['SliceThickness'])

#     # if you prefer to pull the information you need
#     pt = df_cirrus['PC_NUMBER'].iloc[ix] 
#     ndr = df_cirrus['NoduleID'].iloc[ix]
#     px = float(df_cirrus['PixelSpacing'].iloc[ix].split('\\')[-1])
#     sl_thicc = float(df_cirrus['SliceThickness'].iloc[ix])

#     print(f'Working on Patient#{pt} for nodule #{ndr}...')

#     # kernel type:
#     kern = df_cirrus['ConvolutionKernel'].iloc[ix]

#     xc = df_cirrus.CoordX.iloc[ix]
#     yc = df_cirrus.CoordY.iloc[ix]
#     zc = df_cirrus.CoordZ.iloc[ix]

#     # get the minor and major axis length
#     min_axis = df_cirrus.MinorAxis_mm.iloc[ix] / 2 # radius
#     max_axis = df_cirrus.MajorAxis_mm.iloc[ix] / 2 

#     BBOX_Px_Z_CENTER = df_cirrus['SliceNumber'].iloc[ix] + 1 # ASSUMPTION: Cirrus uses 0-indexing, hence the +1
#     BBOX_Z = [np.floor(BBOX_Px_Z_CENTER - (max_axis / px)), np.ceil(BBOX_Px_Z_CENTER + (max_axis / px))]

#     # this method is not efficient as it will reload all the dicoms for each nodule, despite being from the same patient
#     path_to_dicom = rf'..\\path\\to\\{pt}\\dicom\\'
#     dicom_list = walkDir(path_to_dicom)

#     # iterate through dicoms, find what convolution kernel you need
#     volume = 0
#     sl_last = 0
#     for fp in dicom_list:
#         # try to read dicom
#         try:
#             dcm = pyd.dcmread(fp)
#         except pyd.InvalidDicomError as invDError:
#         #   print(invDError)
#             continue

#         # try to read convolution kernel and image type
#         try:
#             conv_kernel = dcm.ConvolutionKernel  # this "filter" is a flawed process. See comment in error.
#             img_type = str(dcm.ImageType)
#         except AttributeError as e:
#             print(e)  # This should not trigger by the time we're done as the BBOX might be need to be laid over CXR too

#         # only keeps images with matching convolution kernel
#         if conv_kernel != kern:
#             continue

#         dcm_im = dcm.pixel_array # should be the dicom image, likley in Hounsfield Units (HU)
#         # it may also be the 'Pixel Data' tag
#         # if you want to visualize the image, you will have to add the RescaleIntercept value to the image
#         # and use the Lung window value, default is [-600, 1500] I believe, and then finally convert it to uint8 [px val 0-255]

#         if not volume:
#             im_shape = np.shape(dcm_im)
#             volume = np.zeros((im_shape[0],im_shape[1],2000))
        
#         dcm_keep = dcm
#         PATIENT_POSITION = dcm['ImagePositionPatient'].value
#         sl = dcm.InstanceNumber - 1 # or sometimes it is dcm.SliceNumber. Can't remember if this is 0-indexed or not, I'm assuming no, hence the -1.
#         volume[:,:,sl] = dcm_im[:,:,:] # this might give you some headaches. You could vectorize it all then reshape back too.
        
#         if sl > sl_last:
#             sl_last = sl # keeping the last slice_value

#         xc_conv = coord_conversion(xc, PATIENT_POSITION[0], px) # gets real world distance between xc and patient position
#         yc_conv = coord_conversion(yc, PATIENT_POSITION[1], px)
#         zc_conv = coord_conversion(zc, PATIENT_POSITION[2], sl_thicc)

#         # I'm guessing here with this
#         if BBOX_Z[0] <= zc_conv <= BBOX_Z[1]:
#             do_something = 1
#             #here we would maybe have to include the Graphics Annotation Tag for this particular slice

#     volume = volume[:,:,0:sl_last] # there is a more elegant way of trimming the excess zeroes off the volume for viewing, but this gets the job done. 
    




# chatgpt suggested sr worflow

# import pydicom
# from datetime import datetime
# from highdicom.sr import ComprehensiveSR, ImageRegion3D
# from highdicom.sr.coding import CodedConcept
# from highdicom.content import GraphicTypeValues, ImageLibraryEntryDescriptor
# import uuid

# # Simulated image metadata (you should replace these with real DICOM values)
# sop_instance_uids = [
#     '1.2.3.4.5.6.7.1',
#     '1.2.3.4.5.6.7.2',
#     '1.2.3.4.5.6.7.3',
# ]
# study_instance_uid = '1.2.3.4.5.6'
# series_instance_uid = '1.2.3.4.5.6.1'
# frame_of_reference_uid = '1.2.3.4.5.6.99'

# # Create image library references
# image_library = [
#     ImageLibraryEntryDescriptor(
#         sop_class_uid='1.2.840.10008.5.1.4.1.1.2',  # CT Image Storage
#         sop_instance_uid=uid,
#         study_instance_uid=study_instance_uid,
#         series_instance_uid=series_instance_uid,
#         frame_of_reference_uid=frame_of_reference_uid
#     )
#     for uid in sop_instance_uids
# ]

# # Example 3D bounding box (polyline of corners, 5 bottom, 5 top)
# graphic_data = [
#     (100.0, 120.0, 50.0),
#     (110.0, 120.0, 50.0),
#     (110.0, 130.0, 50.0),
#     (100.0, 130.0, 50.0),
#     (100.0, 120.0, 50.0),
#     (100.0, 120.0, 60.0),
#     (110.0, 120.0, 60.0),
#     (110.0, 130.0, 60.0),
#     (100.0, 130.0, 60.0),
#     (100.0, 120.0, 60.0),
# ]

# # Create a 3D region object
# region = ImageRegion3D(
#     graphic_type=GraphicTypeValues.POLYLINE,
#     graphic_data=graphic_data,
#     frame_of_reference_uid=frame_of_reference_uid
# )

# # Build the SR document
# sr = ComprehensiveSR(
#     evidence=[],
#     observer_context=[],
#     procedure_reported=CodedConcept('125203', 'DCM', 'Measurement Report'),
#     title=CodedConcept('125203', 'DCM', 'Measurement Report'),
#     series_instance_uid=str(uuid.uuid4()),
#     series_number=1,
#     instance_number=1,
#     sop_instance_uid=str(uuid.uuid4()),
#     manufacturer='MyApp',
#     content_date=datetime.now().date(),
#     content_time=datetime.now().time(),
#     image_library=image_library,
#     content=[region]
# )

# # Save to file
# sr.save_as('example_sr_boundingbox.dcm')





# rest of process: 

# 1. create a test for dicom files to verify the filters of the dicom file loader

# 2. figure out what highdicom.sr.3DSR or whatever it's called needs in its inputs, 
# ESPECIALLY the number of image inputs (if even necessary) it needs to overlay on -
# does it need images to overlay on at all or does it just mark down the bounding box
# in reference to an image? If the latter, does it need to verify whether the bounding
# box is accurate or not?

# 3. try to figure out how patient perspective as well as pixel coordinates vs coordinates
# relative to patient and ensure that the CIRRUS output can be converted into the
# desired outputs before entering these pieces of information into the system
