#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
indep_functions
To do:
    non-object-based independent functions to support hipims IO
-----------------    
Created on Thu Apr 23 11:45:24 2020

@author: Xiaodong Ming
"""
import gzip
import pickle
import os
import shutil
import numpy as np
from .rainfall_processing import _check_rainfall_rate_values

def load_object(file_name):
    """ Read a pickle file as an InputHipims/OutputHipims object
    """
    #read an object file
    try:
        with gzip.open(file_name, 'rb') as input_file:
            obj = pickle.load(input_file)
    except:
        with open(file_name, 'rb') as input_file:
            obj = pickle.load(input_file)
    print(file_name+' loaded')
    return obj

def save_object(obj, file_name, compression=True):
    """ Save an InputHipims/OutputHipims object to a pickle file 
    """
    # Overwrites any existing file.
    if compression:
        with gzip.open(file_name, 'wb') as output_file:
            pickle.dump(obj, output_file, pickle.HIGHEST_PROTOCOL)
    else:
        with open(file_name, 'wb') as output_file:
            pickle.dump(obj, output_file, pickle.HIGHEST_PROTOCOL)
    print(file_name+' has been saved')

def clean_output(case_folder, num_of_sections, file_tag='*'):
    """ delete contents in output folder(s)
    """
    if case_folder[-1]!='/':
            case_folder = case_folder+'/'
    if num_of_sections==1:
        files_to_remove = case_folder+'/output/'+file_tag
        os.system('rm '+files_to_remove)
    else:    
        for i in range(num_of_sections):
            files_to_remove = case_folder+str(i)+'/output/'+file_tag
            os.system('rm '+files_to_remove)

#%% ***************************************************************************
# *************************Public functions************************************
def write_times_setup(case_folder=None, num_of_sections=1, time_values=None):
    """
    Generate a times_setup.dat file. The file contains numbers representing
    the start time, end time, output interval, and backup interval in seconds
    time_values: array or list of int/float, representing time in seconds,
        default values are [0, 3600, 1800, 3600]
    """
    case_folder = _check_case_folder(case_folder)
    if time_values is None:
        time_values = np.array([0, 3600, 1800, 3600])
    time_values = np.array(time_values)
    time_values = time_values.reshape((1, time_values.size))
    if num_of_sections == 1:
        np.savetxt(case_folder+'/input/times_setup.dat', time_values, fmt='%g')
    else:
        np.savetxt(case_folder+'/times_setup.dat', time_values, fmt='%g')
    print('times_setup.dat created')

def write_device_setup(case_folder=None,
                       num_of_sections=1, device_values=None):
    """
    Generate a device_setup.dat file. The file contains numbers representing
    the GPU number for each section
    case_folder: string, the path of model
    num_of_sections: int, the number of GPUs to use
    device_values: array or list of int, representing the GPU number
    """
    case_folder = _check_case_folder(case_folder)
    if device_values is None:
        device_values = np.array(range(num_of_sections))
    device_values = np.array(device_values)
    device_values = device_values.reshape((1, device_values.size))
    if num_of_sections == 1:
        np.savetxt(case_folder+'/input/device_setup.dat',
                   device_values, fmt='%g')
    else:
        np.savetxt(case_folder+'/device_setup.dat', device_values, fmt='%g')
    print('device_setup.dat created')

def write_rain_source(rain_source, case_folder=None, num_of_sections=1):
    """ Write rainfall sources [Independent function from hipims class]
    rain_source: numpy array, The 1st column is time in seconds, the 2nd
        towards the end columns are rainfall rate in m/s for each source ID in
        rainfall mask array
    if for multiple GPU, then copy the rain source file to all domain folders
    case_folder: string, the path of model
    """
    rain_source = np.array(rain_source)
    # check rainfall source value to avoid very large raifall rates
    _ = _check_rainfall_rate_values(rain_source)
    case_folder = _check_case_folder(case_folder)
    fmt1 = '%g'  # for the first col: times in seconds
    fmt2 = '%.8e'  # for the rest array for rainfall rate m/s
    num_mask_cells = rain_source.shape[1]-1
    format_spec = [fmt2]*num_mask_cells
    format_spec.insert(0, fmt1)
    if num_of_sections == 1: # single GPU
        file_name = case_folder+'input/field/precipitation_source_all.dat'
    else: # multi-GPU
        file_name = case_folder+'0/input/field/precipitation_source_all.dat'
    with open(file_name, 'w') as file2write:
        file2write.write("%d\n" % num_mask_cells)
        np.savetxt(file2write, rain_source, fmt=format_spec, delimiter=' ')
    if num_of_sections > 1:
        for i in np.arange(1, num_of_sections):
            field_dir = case_folder+str(i)+'/input/field/'
            shutil.copy2(file_name, field_dir)
    print('precipitation_source_all.dat created')

#%% =======================Value check functions===============================
def _check_case_folder(case_folder):
    """ check the format of case folder
    """
    if case_folder is None:
        case_folder = os.getcwd()
    if not case_folder.endswith('/'):
        case_folder = case_folder+'/'
    return case_folder