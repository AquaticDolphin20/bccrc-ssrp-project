# bccrc-ssrp engine base, by: Ian Janzen, 2025-06-10

import os
import pandas as pd
from pathlib import Path
import numpy as np
import sqlite3 as sql
import time
# from scipy.spatial import distance
import matplotlib.pyplot as plt

import pydicom as pyd 
from pydicom.errors import InvalidDicomError
from pydicom.misc import is_dicom
from pydicom.sr.codedict import codes
import highdicom as hd
from highdicom.sr import CodedConcept
from pydicom.uid import generate_uid

# code is split into sections which contain:
# 1. functions to find dicoms using cirrus outputs, filtering for study, relevant slices of ROI, and reconstruction kernel
# 2. functions to create an ellipsoid annotation using cirrus outputs as a segmentation mask, which is to be inputted into SR file
# 3. functions to create or search for the appropriate uids to input within the SR file, done through maintaining a local sql database
# 4. an overarching loop that for each cirrus row finds dicoms, creates the annotation, finds/creates uids, and outputs this annotation as an SR file


# cirrus_dictionary contains keys as variables within the engine that can be called to return a value, which represents a cirrus column name.
# thus, if the cirrus column outputs change, please change the values to the corresponding new column output name.
# for example, "patient_id", which is a variable used inside this program, is mapped to "PC_NUMBER", which is the name of the cirrus column.
cirrus_dictionary = {
    "patient_id": "PC_NUMBER",
    "nodule_major_axis": "MajorAxis",
    "nodule_minor_axis": "MinorAxis",
    "nodule_diameter": "Diameter_mm"
}

# but what about parsing? think about this later


# reads in all dicoms - we technically don't need this function at all within the workflow - this function is obsolete and for test purposes only
# the reason we don't want to read in all files is to decrease runtime
# def load_dicom_old(dicom_path):
#     dicom_datasets = [] 
#     for filename in os.listdir(dicom_path):
#         filepath = os.path.join(dicom_path, filename)

#         # filter for DICOM files add on: 
#         # if all dicom files end with .dcm, this is a useful filter to decrease runtime
#         if not filename.lower().endswith('.dcm'): 
#             continue

#         # filter for DICOM file option 1
#         try:
#             ds = pyd.dcmread(filepath) # can add argument stop_before_pixels = True to dcmread if we don't need to load all images
#             dicom_datasets.append(ds)
#         except InvalidDicomError as e: # if not dicom file, print error
#             print("Invalid DICOM: ", e)
#             continue

#         # filter for DICOM file option 2
#         # if is_dicom(filepath): 
#             # ds = pyd.dcmread(filepath)

#     return dicom_datasets # returns datasets, not filepaths

# dicom_data = load_dicom_old('../test_data/dicom_data/52100009_baseline')

# def fetch_uids_old(database, study_uid, patient_id, nodule_id)

#     cursor = database.cursor()

#     # fetch tracking_uid if possible
#     cursor.execute('select tracking_uid from uid_list where patient_id = ? and nodule_id = ?', 
#                     (patient_id, nodule_id))
#     tracking_uid = cursor.fetchone()
    
#     # fetch series_uid if possible
#     cursor.execute('select sr_series_uid from uid_list where study_uid = ?',
#                     (study_uid,))
#     sr_series_uid = cursor.fetchone()

#     if tracking_uid and sr_series_uid: # line prevents crashing if there no rows that match criteria
#         return [tracking_uid[0], sr_series_uid[0]]

#     else: 
#         return None



# reads in cirrus csv - created a temporary read in, need to fix parsing
cirrus_data = pd.read_csv('../test_data/cirrus_data/PanCan_SampleDoc_CIRRUS.csv') # tested and works

# path_to_csv = r'..\\path\\to\\csv\\..'

# df_cirrus = pd.read_csv(path_to_csv)



# # functions below are involved with bounding box creator and filtering process

def load_dicoms(dicom_path, dicom_attr: str, desired_values: list[str], need_image: bool = True):
    # helper function for find_slices
    # desired_values is a list just in case we ever need multiple options for filtering, ex want to filter for multiple studies
    # returns dicom images from a database that has desired (can change input) dicom_attr match desired_values, which is a list of strings
    # filters dicoms from a dicom folder path, the filter_dicoms function filters from a list of dicoms

    filtered_dicom_datasets = [] # list will hold dicom folders matching desired dicom_attr
    for root, dirs, files in os.walk(dicom_path):
        for filename in files:
            filepath = os.path.join(root, filename)

            # filter for DICOM files add on: 
            # if all dicom files end with .dcm, this is a useful filter to decrease runtime
            if not filename.lower().endswith('.dcm'): 
                continue

            # filter for DICOM file option 1
            try:
                if need_image == True:
                    ds = pyd.dcmread(filepath)
                else: 
                    ds = pyd.dcmread(filepath, stop_before_pixels = False) # reduce run time by not loading images where image loading isn't necessary

                # filter process
                if getattr(ds, dicom_attr, None) in desired_values: # should we use tags instead of dicom attribute name (more prone to change)
                    filtered_dicom_datasets.append(ds) 

            except InvalidDicomError as e: # if not dicom file, print error
                print("Invalid DICOM: ", e)
                continue

    return filtered_dicom_datasets

def filter_dicoms(dicoms: list, dicom_attr: str, desired_values: list[str]):
    # helper function for find_slices
    # returns dicom images from a database that has desired (can change input) dicom_attr match desired_values, which is a list of strings
    # filters dicoms from a list of dicoms

    filtered_dicom_datasets = []
    for dicom in dicoms:
        if getattr(dicom, dicom_attr, None) in desired_values: # should we use tags.. 
            filtered_dicom_datasets.append(dicom)

    return filtered_dicom_datasets

def find_slice_numbers(df_cirrus, row):
    # helper function for find_slices
    min_slice = int(np.floor(df_cirrus.SliceNumber.iloc[row] - (df_cirrus.Diameter_mm.iloc[row] / 2) / df_cirrus.SliceThickness.iloc[row])) # do we need to convert SliceNumber, Diameter_mm, and SliceThickness into floats
    max_slice = int(np.ceil(df_cirrus.SliceNumber.iloc[row] + (df_cirrus.Diameter_mm.iloc[row] / 2) / df_cirrus.SliceThickness.iloc[row]))
    target_slices = list(range(min_slice, max_slice + 1)) # +1 since range goes up to but not including the number of the max slice 
    
    return target_slices

def find_slices(dicom_path, df_cirrus, row, need_images: bool = True):
    # high level dicom filtering function that 
    # 1. finds the right dicoms for the cirrus patient row
    # 2. finds the right slices for the cirrus patient row
    # 3. finds the right kernel for the cirrus patient row
    # then returns a list of dicoms specific to the cirrus output
    
    # old code
    # patient_identifier = [str(int(df_cirrus.PC_NUMBER.iloc[row]))] # pc number will give you multiple studies of the same patient instead of just one study
    
    # suggested filter to load the study
    patient_identifier = [str(df_cirrus.StudyUID.iloc[row])] # there is currently no cirrus output for study UID, row has been manually updated in test data

    # first filter returns only dicoms with the same patientName or a similar identifier (we can change it)
    patient_dicoms = load_dicoms(dicom_path, "StudyInstanceUID", patient_identifier, need_images)

    # finds slices numbers relevant to annotations 
    patient_slice_numbers = find_slice_numbers(df_cirrus, row) # patient_slice_numbers should hold a list of ints
    
    # second filter returns only target slices from dicom database involved with bbox
    patient_slices = filter_dicoms(patient_dicoms, "InstanceNumber", patient_slice_numbers)

    # third filter returns only target slices with kernel - three separate filters are not ideal in terms of UI or cleanliness but it does make the code easier to change since the functions are less complicated
    # we want annotations to appear in different kernels as well? If so just make patient_slices_kernel a list of list of slices per kernel
    patient_slices_kernel = filter_dicoms(patient_slices, "ConvolutionKernel", df_cirrus.ConvolutionKernel.iloc[row])
  
    return patient_slices_kernel

def create_annotation_point(df_cirrus, row):
    # creates a point annotation based on nodule center output in a cirrus row - this is eligible to put into ImageRegion3D

    x_center = df_cirrus.CoordX.iloc[row]
    y_center = df_cirrus.CoordY.iloc[row]
    z_center = df_cirrus.CoordZ.iloc[row]

    # create numpy object
    annotation_numpy = np.array([(x_center, y_center, z_center)])

    return annotation_numpy

def create_annotation_ellipsoid(df_cirrus, row):
    # creates an ellipsoid based on the information in a cirrus row
    # an ellipsoid is a 3d shape defined by a numpy array, containing six numpy 3d coordinates with each having an (x, y, z) defined in patient space
    # 1. the first set of two coordinates (coordinates 1 and 2) in the array specify extremities for axis x
    # 2. the second set of two coordinates (coordinates 3 and 4) in the array specify extremities for axis y
    # 3. the third set of two coordinates (coordinates 5 and 6) in the array specify extremities for axis z

    # take information for easier access - values should be in patient position already, no need to divide by pixel spacing
    x_center = df_cirrus.CoordX.iloc[row]
    y_center = df_cirrus.CoordY.iloc[row]
    z_center = df_cirrus.CoordZ.iloc[row]
    maj_axis = df_cirrus.MajorAxis_mm.iloc[row] / 2 # divided by two since maj axis is a diameter
    radius = df_cirrus.Diameter_mm.iloc[row] / 2 # original code did not divide diameter by 2

    # compute ellipsoid coordinate edges
    x_min, x_max = np.floor(x_center - maj_axis), np.ceil(x_center + maj_axis)
    y_min, y_max = np.floor(y_center - maj_axis), np.ceil(y_center + maj_axis)
    z_min, z_max = z_center - radius, z_center + radius # do these rows need np.floor or np.ceil to round to be safe or is cirrus diameter accurate
    
    # create ellipsoid numpy object using coordinate edges and center
    annotation_numpy = np.array([
        (x_min, y_center, z_center), (x_max, y_center, z_center), # extremities on x axis
        (x_center, y_min, z_center), (x_center, y_max, z_center), # extremities on y axis
        (x_center, y_center, z_min), (x_center, y_center, z_max) # extremities on z axis
    ])

    return annotation_numpy


# ---------------------------------------------------------------------------------------


def find_pixel_coordinates(world_coord: np.ndarray, dicom_slices, df_cirrus, row):
    # takes center from cirrus output -> converts to pixel position

    # get information ready for conversion
    image_orientation_patient = getattr(dicom_slices[0], "ImageOrientationPatient", None) # example: [1, -7.04487467e-016, 0, 0, 0, -1], only need this for non-axial scans (eg oblique scans)
    pixel_spacing_x, pixel_spacing_y = map(float, df_cirrus.PixelSpacing.iloc[row].split('\\'))
    slice_thickness = float(df_cirrus.SliceThickness.iloc[row]) # example: 1

    # get patient image position from dicom
    image_position = np.array(getattr(dicom_slices[0], "ImagePositionPatient", None)) # example: [-255.5, -132, -1.5]

    # convert for each axis
    voxel_coord = np.array([
        (world_coord[0] - image_position[0]) / pixel_spacing_x, # x voxel
        (world_coord[1] - image_position[1]) / pixel_spacing_y, # y voxel
        (image_position[2] - world_coord[2]) / slice_thickness # z voxel - because feet first scan, slices are feet to head
        ])

    return voxel_coord

def create_segmentation_ellipsoid(dicom_slices, df_cirrus, row):
    
    # get world coordinates np array ellipsoid and point (center)
    annotation_ellipsoid = create_annotation_ellipsoid(df_cirrus, row)
    annotation_point = create_annotation_point(df_cirrus, row)

    # turn each point of ellipsoid into pixel coordinates
    annotation_ellipsoid_pixel = []

    for point in annotation_ellipsoid:
        pixel_ellipsoid_point = find_pixel_coordinates(point, dicom_slices, df_cirrus, row)
        annotation_ellipsoid_pixel.append(pixel_ellipsoid_point)

    annotation_ellipsoid_np = np.array(annotation_ellipsoid_pixel)

    # turn each point of point into pixel coordinates
    annotation_center = find_pixel_coordinates(annotation_point[0], dicom_slices, df_cirrus, row)

    # create ellipsoid segmentation binary mask
    volume_shape = (len(dicom_slices), dicom_slices[0].Rows, dicom_slices[0].Columns) # z, y, x volume shape

    maj_axis = df_cirrus.MajorAxis_mm.iloc[row] / 2 # divided by two since maj axis is a diameter
    
    # clean up code below put it in another function
    pixel_spacing_x, pixel_spacing_y = map(float, df_cirrus.PixelSpacing.iloc[row].split('\\'))
    slice_thickness = float(df_cirrus.SliceThickness.iloc[row])

    # compute radii
    radii_x = (df_cirrus.MajorAxis_mm.iloc[row] / 2) / pixel_spacing_x
    radii_y = (df_cirrus.MajorAxis_mm.iloc[row] / 2) / pixel_spacing_y
    radii_z = (df_cirrus.Diameter_mm.iloc[row] / 2) / slice_thickness
    
    z, y, x = np.ogrid[0:volume_shape[0], 0:volume_shape[1], 0:volume_shape[2]] # volume shape is z, y, x
    # z becomes array([0, 1, 2, ... 510, 511]) if the volume_shape was 512
    # these three coordinates should cover the entire image as a mask, so the mask is huge (covers slice) but ellipsoid only covers tiny bit of mask
    
    segmentation_mask = (
        ((z - annotation_center[2]) / radii_z) ** 2 + # essentially writing ((512 - 210) / 3)**2
        ((y - annotation_center[1]) / radii_y) ** 2 +
        ((x - annotation_center[0]) / radii_x) ** 2
    ) <= 1 # each point on the grid will have its (x, y, z) plugged into the mathematical equation above, and if the equation computes to <= 1 it will return False
    # thus each point within the ellipsoid will <= 1 and each point outside the ellipsoid will be > 1.

    return segmentation_mask



# sql database initializer 
def load_uids(database_path = "../sql_output/nodule_uids.db"): # change default path for parsing later
    
    # forms connection object with database in path
    conn = sql.connect(database_path)
    cursor = conn.cursor()
    
    # creates table called patient_uids if it doesn't already exist in the db file (will do this first time running python script)
    cursor.execute('''
        create table if not exists uid_list ( -- creates a new table called "uids" if it doesn't already exist
            study_uid text not null, 
            study_date text not null,
            patient_id text not null, 
            nodule_id text not null, 
            tracking_uid text not null, 
            seg_series_uid text not null,
            seg_sop_uid text not null,
            sr_series_uid text not null, 
            sr_sop_uid text not null,
            PRIMARY KEY (study_uid, patient_id, nodule_id) -- prevents duplicate entries for same patient and nodule id
        )
    ''')
    conn.commit()

    return conn

# tracking_uid adder for new tracking uids
def save_uids(database, study_uid: str, study_date: str, patient_id: str, nodule_id: str, tracking_uid: str, seg_series_uid: str, seg_sop_uid: str, sr_series_uid: str, sr_sop_uid: str):
    cursor = database.cursor()
    cursor.execute('insert or ignore into uid_list (study_uid, study_date, patient_id, nodule_id, tracking_uid, seg_series_uid, seg_sop_uid, sr_series_uid, sr_sop_uid) values (?, ?, ?, ?, ?, ?, ?, ?, ?)',
                    (study_uid, study_date, patient_id, nodule_id, tracking_uid, seg_series_uid, seg_sop_uid, sr_series_uid, sr_sop_uid))
    database.commit()

# uid_fetcher helper function
def fetch_uid(database, fetch_columns: str, filter_columns: list[str], filter_values: list):

    cursor = database.cursor()

    # build sql query 
    where_clause = ' and '.join([f"{column} = ?" for column in filter_columns])
    query = f"select {fetch_columns} from uid_list where {where_clause}"

    # execute query
    cursor.execute(query, tuple(filter_values))
    uid = cursor.fetchone()

    if uid: # line prevents errors if there no rows that match criteria, ie new uid must be created
        return uid[0]

    else: 
        return None


# top level function that saves everything if it isn't already present and returns dictionary with nodule id, tracking UID for a specific cirrus row and column
def find_uids(df_cirrus, row): 
    # helper functions are initialize_database, loading_tracking_uids_as_dict, etc

    # load database
    database = load_uids() # initializes database
    
    # get information from cirrus
    study_uid = str(df_cirrus.StudyUID.iloc[row])
    study_date = str(df_cirrus.StudyDate.iloc[row])
    patient_id = str(int(df_cirrus.PC_NUMBER.iloc[row]))
    nodule_id = str(int(df_cirrus.LesionID.iloc[row])).zfill(4)

    # fetch uids if possible
    tracking_uid = fetch_uid(database, "tracking_uid", ["patient_id", "nodule_id"], [patient_id, nodule_id]) # tracking_id is based on patient, nodule, same across different studies over time
    seg_series_uid = fetch_uid(database, "seg_series_uid", ["study_uid"], [study_uid])
    seg_sop_uid = fetch_uid(database, "seg_sop_uid", ["study_uid", "patient_id", "nodule_id"], [study_uid, patient_id, nodule_id])
    sr_series_uid = fetch_uid(database, "sr_series_uid", ["study_uid"], [study_uid]) # series_uid is based on patient, date, different across different studies over time
    sr_sop_uid = fetch_uid(database, "sr_sop_uid", ["study_uid", "patient_id", "nodule_id"], [study_uid, patient_id, nodule_id])

    # create tracking uid if nonexistant for a nodule + patient
    if not tracking_uid: # uid will return True on if statement if uid is NoneType ie nonexistant
        tracking_uid = generate_uid(prefix = '1.2.826.0.1.3680043.8.498.') # prefix is DICOM test prefix, please change to BC Cancer organization prefix

    # segmentation uid creation if necessary
    # i am aware that sometimes the program you are using to segment files already has a series instance uid/sop instance uid, so just input that accordingly through changing the code
    if not sr_sop_uid:
        seg_series_uid = generate_uid(prefix = '1.2.826.0.1.3680043.8.498.') # prefix is DICOM test prefix, please change to segmentator program prefix
    
    if not sr_sop_uid:
        seg_sop_uid = generate_uid(prefix = '1.2.826.0.1.3680043.8.498.') # prefix is DICOM test prefix, please change to segmentator program prefix
    
    
    # sr uid creation if necessary
    # create sr_series_uid if nonexistant for a study
    if not sr_series_uid: 
        sr_series_uid = generate_uid(prefix = '1.2.826.0.1.3680043.8.498.') # prefix is DICOM test prefix, please change to BC Cancer organization prefix
    
    # create sr_instance_uid - new instance is created every time BUT
    # we check first just in case this program was run accidentally already so we get the sop_uid that's already written in db
    if not sr_sop_uid:
        sr_sop_uid = generate_uid(prefix = '1.2.826.0.1.3680043.8.498.') # prefix is DICOM test prefix, please change to BC Cancer organization prefix
    
    # save uids to database
    save_uids(database, study_uid, study_date, patient_id, nodule_id, tracking_uid, seg_series_uid, seg_sop_uid, sr_series_uid, sr_sop_uid) # registers unique uids into sql file for future use

    database.close()

    return [tracking_uid, seg_series_uid, seg_sop_uid, sr_series_uid, sr_sop_uid]


# cirrus loop where each row will create the bounding box for that specific row of information
# for row in range(len(cirrus_data))[:1]: # the [:1] is only for tests

for row in [30]:
    # get nodule number for file naming (ex. sr_1.dcm, sr_2.dcm depending on nodule)
    nodule_number = str(int(cirrus_data.LesionID.iloc[row])).zfill(4) # zfill allows at least four zeros in front of nodule number, so "0004" instead of "4"

    # get and save tracking UID
    tracking_uid = find_uids(cirrus_data, row)[0]
    seg_series_uid = find_uids(cirrus_data, row)[1]
    seg_sop_uid = find_uids(cirrus_data, row)[2]
    sr_series_uid = find_uids(cirrus_data, row)[3]
    sr_sop_uid = find_uids(cirrus_data, row)[4]
    

    placeholder_uid = "1.2.826.0.1.3680043.8.498.183683720340460575189739599081"

    # get list of dicom slices involved in ROI - need to create parsing for filepath
    slices = find_slices('../test_data/dicom_data', cirrus_data, row)

    # get a list of SourceImageForSegmentation objects
    source_images = []
    
    for slice in slices: # change variable name to something other than "slice"
        image = hd.sr.SourceImageForSegmentation(
            referenced_sop_class_uid = getattr(slice, 'SOPClassUID', None),
            referenced_sop_instance_uid = getattr(slice, 'SOPInstanceUID', None)
        )

        source_images.append

    # get numpy annotation
    annotation_point = create_annotation_point(cirrus_data, row)
    annotation_ellipsoid = create_annotation_ellipsoid(cirrus_data, row)
    segmentation_ellipsoid = create_segmentation_ellipsoid(slices, cirrus_data, row)

    # TESTS for slices
    # print(len(slices))
    # print(getattr(slices[0], 'PatientName', None))
    # print(getattr(slices[-1], 'PatientName', None))
    # print(getattr(slices[0], 'InstanceNumber', None))
    # print(getattr(slices[-1], 'InstanceNumber', None))
    # print(getattr(slices[0], 'ConvolutionKernel', None))
    # print(getattr(slices[-1], 'ConvolutionKernel', None))
    # for ds in slices: 
    #     print(ds.filename)

    # TESTS for annotation
    # print(annotation_point)
    # print(annotation_ellipsoid)
    # print(segmentation_ellipsoid)

    # TESTS for uids
    # print(tracking_uid)
    # print(sr_series_uid)
    # print(sr_sop_uid)


# ------------------------------------------------------------------------------------



    # maybe don't need this: 
    region_point = hd.sr.ImageRegion3D(
        graphic_type = hd.sr.GraphicTypeValues3D.POINT, 
        graphic_data = annotation_point,
        frame_of_reference_uid = getattr(slices[0], "FrameOfReferenceUID", None)
    )

    # print(region_point)

    region_ellipsoid = hd.sr.Scoord3DContentItem(
        name = CodedConcept(value = '111030', # value for ROI
                            scheme_designator = 'DCM', 
                            meaning = 'Region of Interest'), # CodedConcept must be manually inputted; existing patient dicoms do not have such a field
        graphic_type = hd.sr.GraphicTypeValues3D.ELLIPSOID,
        graphic_data = annotation_ellipsoid, 
        frame_of_reference_uid = getattr(slices[0], "FrameOfReferenceUID", None), # can pull frameofreferenceUID from first slice since it doesn't change per scan
    )

    # print(region_ellipsoid)

    segment_description = hd.seg.SegmentDescription(
        segment_number = 1, # replace later
        segment_label = "Nodule",
        segmented_property_category = CodedConcept(
            value = 'M-01000',
            scheme_designator = 'SRT',
            meaning = 'Morphologically Abnormal Structure'), # check fields for accuracy later and compare with: https://dicom.nema.org/medical/dicom/current/output/chtml/part16/sect_CID_7150.html
        segmented_property_type = CodedConcept(
            value = 'M-03010',
            scheme_designator = 'SRT',
            meaning = 'Nodule'), # check fields for accuracy later and compare with: https://dicom.nema.org/medical/dicom/current/output/chtml/part16/sect_CID_7151.html
        algorithm_type = "MANUAL", # change to AUTOMATIC later and do the algorithm_identification
        # algorithm_identification = "EllipsoidConverter" # optional but we should probably put something here
    )

    # print(segment_description)

    segment = hd.seg.Segmentation(
        source_images = slices,
        pixel_array = segmentation_ellipsoid,
        segmentation_type = "BINARY",
        segment_descriptions = [segment_description], # change later
        series_instance_uid = seg_series_uid, # change later
        series_number = 1, # change later
        sop_instance_uid = seg_sop_uid, # change later
        instance_number = int(nodule_number), # change later
        manufacturer = "Cirrus", # change later
        manufacturer_model_name = "Cirrus", # change later
        software_versions = ("Cirrus v1.0.0", "Cirrus v1.0"), # change later
        device_serial_number = "CIRRUS-v1.0.0-20240101", # change later
        content_description = "Test Segmentation"
    )

    # print(segment)

    segment.save_as(f"../test_output/test8/seg_{int(nodule_number)}.dcm")


    volumetric_roi_tracking_identifier = hd.sr.TrackingIdentifier(
        identifier = f'NoduleID_{nodule_number}',
        uid = tracking_uid
    )


# --------------------------------------- for ellipsoid vector

#     # print(volumetric_roi_tracking_id)

#     vol_surface = hd.sr.VolumeSurface(
#         graphic_type = hd.sr.GraphicTypeValues3D.ELLIPSOID, 
#         graphic_data = annotation_ellipsoid,
#         frame_of_reference_uid = getattr(slices[0], "FrameOfReferenceUID", None), # can pull frameofreferenceUID from first slice since it doesn't change per scan
#         source_images = source_images
#     )

#     # print(vol_surface)

#     vol_group = hd.sr.VolumetricROIMeasurementsAndQualitativeEvaluations(
#         referenced_volume_surface = vol_surface, 
#         tracking_identifier = volumetric_roi_tracking_identifier
#     )

#     print(vol_group)

# # ----------------------------------------------

#     vol_group = hd.sr.VolumetricROIMeasurementsAndQualitativeEvaluations(
#         referenced_segment = volume_segment,
#         tracking_identifier = volumetric_roi_tracking_identifier
#     )

#     # print(vol_group)

#     referenced_segment = hd.sr.ReferencedSegment(
#         sop_class_uid = "1.2.840.10008.5.1.4.1.1.66.4", # code for segmentation storage
#         sop_instance_uid = seg_sop_uid,
#         segment_number = int(nodule_number)
#     )

#     print(referenced_segment)

#     observer = hd.sr.ObserverContext(
#         observer_type = codes.DCM.Device,
#         observer_identifying_attributes = hd.sr.DeviceObserverIdentifyingAttributes(
#             uid = generate_uid(), # need the cirrus uid from bc cancer
#             manufacturer_name = "Cirrus"
#         )
#     )

#     # print(observer)

#     observation_context = hd.sr.ObservationContext(
#         observer_device_context = observer
#     )

#     # print(observation_context)

#     measurement_report = hd.sr.MeasurementReport(
#         observation_context = observation_context, # not optional need to fill
#         procedure_reported = codes.LN.CTUnspecifiedBodyRegion, # for options, visit: https://dicom.nema.org/medical/dicom/current/output/chtml/part16/sect_CID_100.html
#         imaging_measurements = [vol_group],
#         title = codes.DCM.ImagingMeasurementReport # for options, visit: https://dicom.nema.org/medical/dicom/current/output/chtml/part16/sect_CID_7021.html
#     )

    # print(measurement_report)

    # sr_dataset = hd.sr.Comprehensive3DSR(
    #     evidence = slices, # all datasets referenced in the report
    #     content = measurement_report, 
    #     series_number = 1, # switch later if there are cases where there are multiple sr series from the same study
    #     series_instance_uid = sr_series_uid,
    #     sop_instance_uid = sr_sop_uid, 
    #     instance_number = int(nodule_number), # what do I put for this lol is this supposed to be different for every bounding box per patient
    #     series_description = 'Cirrus Structured Report' 
    # )

    # print(sr_dataset)

    # # sr_dataset.save_as(f"../test_output/test8/sr_{int(nodule_number)}.dcm") # need to fix the path


# -------------------------------------------------------------------------------------


# Ian code

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
    


# --------------------------------------------------------------------------------

# tests

# tests for csv read in
# print(cirrus_data.columns) # ensure all columns are present
# print(cirrus_data.head(3)) # ensure data is organized 

# tests for dicom read in
# for ds in dicom_data[:3]:
#     print(ds.filename) # ensure dicom data is present

# tests for dicom filtering
# test1a = load_dicoms('../test_data/dicom_data/52100009_baseline', "InstanceNumber", [3, 4, 5, 6])
# for ds in test1a[:3]:
#     print(ds.filename)
# print(len(test1a)) # checks number of dicom files in test1a

# test1b = load_dicoms('../test_data/dicom_data', "InstanceNumber", [3])
# for ds in test1b:
#     print(ds.filename)
# print(len(test1b)) # checks number of dicom files in test1b

# test1c = load_dicoms('../test_data/dicom_data', "PatientName", [str(52100660)]) # PatientName is a string
# print(len(test1c))
# print(getattr(test1c[0], 'PatientName', None))
# for ds in test1c: 
#     print(ds.filename)


# tests for dicom filter - test works!
# test2 = load_dicoms('../test_data/dicom_data/52100009_baseline', "InstanceNumber", [3, 4, 5, 6])

# for ds in test2:
#     print(ds.filename)

# test2a = filter_dicoms(test2, "InstanceNumber", [3, 4])
# for ds in test2a:
#     print(ds.filename)

# test2b = filter_dicoms(test2, "SeriesInstanceUID", ["1.3.12.2.1107.5.1.4.50138.30000009041917002767100000000"])
# for ds in test2b:
#     print(ds.filename)

# test2c = filter_dicoms(test2, "ConvolutionKernel", ["B35f"])
# for ds in test2c: 
#     print(ds.filename)


# test to see whether need_image parameter works - it works
# test3a = filter_dicoms('../test_data/dicom_data/52100009_baseline', "InstanceNumber", [3, 4, 5, 6], True)
# test3b = filter_dicoms('../test_data/dicom_data/52100009_baseline', "InstanceNumber", [3, 4, 5, 6], False)

# for i, ds in enumerate(test3a):
#     ds.save_as(f"../test_output/test3a/test3a_{i}.dcm")

# for i, ds in enumerate(test3b):
#     ds.save_as(f"../test_output/test3b/test3b_{i}.dcm")

# not the best designed test but it gets the job done - you can see that if you use need_image = True it includes image and thus
# has a bigger file size when saved - the reason why we can't use 
    # for ds in test3a[:2]:
        # print(os.path.getsize(ds.filename))
# is that it references the original file size, not what is actually in the filtered variable



# test to see whether find_slice_numbers function works - it works
# test4_cirrus_data = pd.read_csv('../test_data/cirrus_data/PanCan_SampleDoc_CIRRUS_test4.csv') # tested and works
# print(cirrus_data.columns) # ensure all columns are present
# print(cirrus_data.head(3))

# print(find_slice_numbers(test4_cirrus_data, 0))
# print(find_slice_numbers(test4_cirrus_data, 1))
# print(find_slice_numbers(test4_cirrus_data, 2))



# test to see whether find_slices function works - WORKS!

# test5_cirrus_data = pd.read_csv('../test_data/cirrus_data/PanCan_SampleDoc_CIRRUS.csv')

# test5a = find_slices('../test_data/dicom_data', test5_cirrus_data, 0)
# print(len(test5a))
# print(getattr(test5a[0], 'PatientName', None))
# print(getattr(test5a[-1], 'PatientName', None))
# print(getattr(test5a[0], 'InstanceNumber', None))
# print(getattr(test5a[-1], 'InstanceNumber', None))
# print(getattr(test5a[0], 'ConvolutionKernel', None))
# print(getattr(test5a[-1], 'ConvolutionKernel', None))
# for ds in test5a: 
#     print(ds.filename)

# test5b = find_slices('../test_data/dicom_data', test5_cirrus_data, 15)
# print(len(test5b))
# print(getattr(test5b[0], 'PatientName', None))
# print(getattr(test5b[-1], 'PatientName', None))
# print(getattr(test5b[0], 'InstanceNumber', None))
# print(getattr(test5b[-1], 'InstanceNumber', None))
# print(getattr(test5b[0], 'ConvolutionKernel', None))
# print(getattr(test5b[-1], 'ConvolutionKernel', None))
# for ds in test5b: 
#     print(ds.filename)
# test 5b DOES NOT WORK since we only have patient 52100048 in dicoms and patient 52100044 in cirrus LOLL

# test5c = find_slices('../test_data/dicom_data', test5_cirrus_data, 30)
# print(len(test5c))
# print(getattr(test5c[0], 'PatientName', None))
# print(getattr(test5c[-1], 'PatientName', None))
# print(getattr(test5c[0], 'InstanceNumber', None))
# print(getattr(test5c[-1], 'InstanceNumber', None))
# print(getattr(test5c[0], 'ConvolutionKernel', None))
# print(getattr(test5c[-1], 'ConvolutionKernel', None))
# for ds in test5c: 
#     print(ds.filename)




# tests to see whether create_annotation_ellipsoid and create_annotation_point function works - WORKS!

# test6_cirrus_data = pd.read_csv('../test_data/cirrus_data/PanCan_SampleDoc_CIRRUS.csv')

# test6a = create_annotation_ellipsoid(test6_cirrus_data, 0)
# print(test6a)

# test6b = create_annotation_ellipsoid(test6_cirrus_data, 30)
# print(test6b) 

# test6c = create_annotation_point(test6_cirrus_data, 0)
# print(test6c)

# test6d = create_annotation_point(test6_cirrus_data, 30)
# print(test6d)




# tests to see whether load_uids, fetch_uid, save_uid, find_uid functions work - Works!

# test7 = load_uids("../sql_output/nodule_uids.db") 
# test7_cirrus = pd.read_csv('../test_data/cirrus_data/PanCan_SampleDoc_CIRRUS.csv')

# test7a = save_uids(test7, "1.2.124.113532.142.71.190.136.20090416.103812.15148805",
#             "20090419", "52100009", "0001", 
#             "1.2.826.0.1.3680043.8.498.73188677548692134182519393705869338291", # tracking_uid
#             "1.2.826.0.1.3680043.8.498.46232019650193572610247613322035420344", # seg series uid
#             "1.2.826.0.1.3680043.8.498.90280968034516665733260840226747440580", # seg sop uid
#             "1.2.826.0.1.3680043.8.498.54740032412729907813776555872269400503", # sr series uid
#             "1.2.826.0.1.3680043.8.498.18368372034046057518973959908118183257" # sr sop uid
#             )

# test7b = save_uids(test7, "1.2.124.113532.142.71.190.136.20090416.103812.15148805",
#             "20090419", "52100009", "0001", 
#             "1.2.826.0.1.3680043.8.498.73188677548692134182519393705869338291", # tracking_uid
#             "1.2.826.0.1.3680043.8.498.46232019650193572610247613322035420344", # seg series uid
#             "1.2.826.0.1.3680043.8.498.90280968034516665733260840226747440580", # seg sop uid
#             "1.2.826.0.1.3680043.8.498.54740032412729907813776555872269400503", # sr series uid
#             "1.2.826.0.1.3680043.8.498.18368372034046057518973959908118183257" # sr sop uid
#             ) # test 7b should not save into nodule_uid, should preserve old uids in test7a due to duplicate keys

# test7c = fetch_uid(test7, "seg_sop_uid", ["study_uid"], ["1.2.124.113532.142.71.190.136.20090416.103812.15148805"])
# print(test7c)

# test7d = fetch_uid(test7, "sr_series_uid", ["study_uid"], ["1.2.124.113532.142.71.190.136.20090416.103812.15148805"])
# print(test7d)

# test7e = fetch_uid(test7, "tracking_uid", ["patient_id", "nodule_id"], ["52100660", "0006"])
# print(test7e) 

# test7f = fetch_uid(test7, "tracking_uid", ["patient_id", "nodule_id"], ["52100300", "0001"])
# print(test7f) # should print None instead of raising error

# test7g = find_uids(test7_cirrus, 0)
# print(test7g)

# test7h = find_uids(test7_cirrus, 30)
# print(test7h)

# test7i = find_uids(test7_cirrus, 36)
# print(test7i)


# find_pixel_coordinates tests - Works!
# test8_cirrus = pd.read_csv('../test_data/cirrus_data/PanCan_SampleDoc_CIRRUS.csv')
# test8_dicom = find_slices('../test_data/dicom_data', test8_cirrus, 0, False)
# test8_world_coord = np.array([
#     test8_cirrus.CoordX.iloc[0],
#     test8_cirrus.CoordY.iloc[0],
#     test8_cirrus.CoordZ.iloc[0]
# ])

# test8a = find_pixel_coordinates(test8_world_coord, test8_dicom, test8_cirrus, 0)
# print(test8a)
# print(test8a[0])
# print(test8a[1])
# print(test8a[2])

# test8b = create_segmentation_ellipsoid(test8_dicom, test8_cirrus, 0)


def visualize_segmentation_mask(mask, num_slices=5):
    # function to visualize segmentation mask for sanity checks

    z_dim = mask.shape[0]
    mid_z = z_dim // 2
    half_window = num_slices // 2

    slice_range = range(max(0, mid_z - half_window), min(z_dim, mid_z + half_window + 1))

    fig, axes = plt.subplots(1, len(slice_range), figsize=(15, 5))
    if len(slice_range) == 1:
        axes = [axes]  # Ensure iterable

    for ax, z_idx in zip(axes, slice_range):
        ax.imshow(mask[z_idx], cmap='gray')
        ax.set_title(f"Z Slice {z_idx}")
        ax.axis('off')

    plt.tight_layout()
    plt.show()

# visualize_segmentation_mask(test8b, num_slices = 2)

# runtime test
# start = time.time()
# conn = load_uids()
# print("Load time:", time.time() - start)


