#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
init
To do:
    initialize a package
Created on Wed Apr  1 14:56:15 2020

@author: Xiaodong Ming
"""
import pkg_resources
import numpy as np
from .InputHipims import InputHipims
from .Raster import Raster

def demo_input(set_example_inputs=True):
    """ A demonstration to generate a hipims input object
    """
    dem_file = pkg_resources.resource_filename(__name__,
                                             'sample/Example_DEM.asc')
    obj_in = InputHipims(dem_data=dem_file)
    if set_example_inputs:
        __set_defaul_input(obj_in)
    # show model summary print(obj_in)
    obj_in.Summary.display()
    fig, ax = obj_in.domain_show(relocate=True, scale_ratio=1000)
    ax.set_title('The Upper Lee catchment')
    return obj_in

def demo_raster(figname=None):
    """ A demonstration to read and show raster files
    figname: the file name to save the figure
    """
    dem_file = pkg_resources.resource_filename(__name__,
                                             'sample/Example_DEM.asc')
    obj_ras = Raster(dem_file)
    fig, ax = obj_ras.mapshow(figname=figname, relocate=True, scale_ratio=1000)
    ax.set_title('The Upper Lee catchment DEM (mAOD)')
    return obj_ras

def get_sample_data():
    """ Get sample data for demonstartion
    Return:
        obj_ras: a DEM raster object
        demo_data: a dictionary with boundary_condition, rain_source, and 
            gauges_pos data
    """
    dem_file = pkg_resources.resource_filename(__name__,
                                             'sample/Example_DEM.asc')
    obj_ras = Raster(dem_file)
    demo_data_file = pkg_resources.resource_filename(__name__,
                                             'sample/Example_data.npy')
    demo_data = np.load(demo_data_file, allow_pickle='TRUE').item()
    return obj_ras, demo_data
    

def __set_defaul_input(obj_in):
    """Set some default values for an InputHipims object
    """
    # load data for the demo
    demo_data_file = pkg_resources.resource_filename(__name__,
                                             'sample/Example_data.npy')
    demo_data = np.load(demo_data_file, allow_pickle='TRUE').item()
    # define initial condition
    h0 = obj_in.Raster.array+0
    h0[np.isnan(h0)] = 0
    h0[h0 < 50] = 0
    h0[h0 >= 50] = 1
    # set initial water depth (h0) and velocity (hU0x, hU0y)
    obj_in.set_parameter('h0', h0)
    obj_in.set_parameter('hU0x', h0*0.0001)
    obj_in.set_parameter('hU0y', h0*0.0002)
    # define boundary condition
    bound_list = demo_data['boundary_condition']
    obj_in.set_boundary_condition(bound_list, outline_boundary='fall')
    # define and set rainfall mask and source (two rainfall sources)
    rain_source = demo_data['rain_source']
    obj_in.set_rainfall(rain_mask=0, rain_source=rain_source)
    # define and set monitor positions
    gauges_pos = demo_data['gauges_pos']
    obj_in.set_gauges_position(gauges_pos)