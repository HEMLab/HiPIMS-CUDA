#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
init
To do:
    initialize a package
Created on Wed Apr  1 14:56:15 2020

@author: Xiaodong Ming
"""
from .demo_functions import demo_input
from .demo_functions import demo_output
from .demo_functions import demo_raster
from .demo_functions import get_sample_data
from .InputHipims import copy_input_obj, InputHipims
from .OutputHipims import OutputHipims
from .indep_functions import load_object, save_object, clean_output
from .indep_functions import write_times_setup, write_device_setup, write_rain_source
