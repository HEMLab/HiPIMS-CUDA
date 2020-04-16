#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
model_summary
To do:
    create, renew, and display model information
Created on Tue Mar 31 16:07:09 2020

@author: Xiaodong Ming
"""
import warnings
import numpy as np
#%% HiPIMS model input information summary
class ModelSummary:
    """ Store and disply all model information including:
    case_folder
    Domain area
    Grid size
    Number of Section
    Initial condition
    Boundary condition
    Rainfall data
    Parameters
    """
    #%======================== initialization function ========================
    def __init__(self, hipims_obj):
        num_valid_cells = hipims_obj._valid_cell_subs[0].size
        self.__case_folder = hipims_obj.case_folder
        num_of_sections = hipims_obj.num_of_sections
        dem_header = hipims_obj.Raster.header
        self.domain_area = num_valid_cells*(dem_header['cellsize']**2)
        runtime_str = hipims_obj.set_runtime(hipims_obj.times)
        birthday_str = hipims_obj.birthday.strftime('%Y-%m-%d %H:%M:%S')
        self.information_dict = {
            '*******************':'Model summary*****************',
            'Case folder':hipims_obj.case_folder,
            'Number of Sections':str(num_of_sections),
            'Grid size':'{:d} rows * {:d} cols, {:g} m resolution'.format(
                dem_header['nrows'], dem_header['ncols'],
                dem_header['cellsize']),
            'Domain area':'{1:,} m^2 with {0:,} valid cells'.format(
                num_valid_cells, self.domain_area),
            'Birthday':birthday_str,
            'Run time':runtime_str,
            '------------':'----Initial Condition----------',
            'h0':'0 for all cells',
            'hU0x':'0 for all cells',
            'hU0y':'0 for all cells',
            'precipitation':'0 for all cells',
            '--------------':'--Boundary Condition---------',
            'Number of boundaries':'1',
            'Boundary details': '(outline) fall, h and hU fixed as zero',
            '---------------':'-Precipitation--------------',
            'precipitation_mask':'0 for all cells',
            'precipitation_source': '0',
            '----------------':'Parameters-----------------'}
        if hasattr(hipims_obj, 'section_id'): # sub-domain object
            self.information_dict['*******************'] = \
                'Section summary****************'
            self.add_items('Domain ID', '{:d}'.format(hipims_obj.section_id))
        else:
            # add all param information if hipims_obj is not a child object
            for key, value in hipims_obj.attributes.items():
                if key not in ['h0', 'hU0x', 'hU0y', 'precipitation',
                               'precipitation_mask', 'precipitation_source']:
                    self.add_param_infor(key, value)

    def display(self):
        """
        Display the model summary information
        """
        for key in self.information_dict.keys():
            if key.endswith('-') or key.endswith('*'):
                print(key+self.information_dict[key])
            else:
                print(key+': '+self.information_dict[key])
        print('*************************************************')

    def write_readme(self, filename=None):
        """
        Write readme file for the summary information
        """
        if filename is None:
            filename = self.__case_folder+'/readme.txt'
        with open(filename, 'w') as file2write:
            for key, value in self.information_dict.items():
                if '------------' in key:
                    file2write.write(key+value+'\n')
                else:
                    file2write.write(key+': '+value+'\n')

    def set_param(self, param_name, param_value):
        """ Set values in the information_dict
        """
        if param_name in self.information_dict.keys():
            self.information_dict[param_name] = param_value
        else:
            raise ValueError(param_name+ ' is not defined in the object')

    def add_items(self, item_name, item_value):
        """ Add a new item to the information_dict
        """
        if type(item_value) is not str:
            item_value = str(item_value)
        self.information_dict[item_name] = item_value

    def add_param_infor(self, param_name, param_value):
        """ Add information of hipims model parameters to the information_dict
        """
        param_value = np.array(param_value)
        item_name = param_name
        if param_value.size == 1:
            item_value = ' {:} for all cells'.format(param_value)
        else:
            if param_name in ['h0', 'hU0x', 'hU0y']:
                num_wet_cells = np.sum(param_value > 0)
                num_wet_cells_rate = num_wet_cells/param_value.size
                item_value = ' nonzero ratio: {:.2f}%'.format(
                    num_wet_cells_rate*100)
            elif param_name == 'precipitation_mask':
                if param_value.dtype != 'int32':
                    param_value = param_value.astype('int32')
                item_value = '{:d} rainfall sources'.format(
                    param_value.max()+1)
            elif param_name == 'precipitation_source':
                rain_max, rain_mean = _check_rainfall_rate_values(param_value)
                display_str = 'max rain: {:.2f} mm/h, '+\
                                'average rain: {:.2f} mm/h'
                item_value = display_str.format(rain_max, rain_mean)
            elif param_name == 'gauges_pos':
                item_value = '{:d} gauges'.format(param_value.shape[0])
            else:
                item_numbers, item_number_counts = \
                    np.unique(param_value, return_counts=True)
                ratio = item_number_counts*100/param_value.size
                formatter = {'float_kind':lambda ratio: "%.3f" % ratio}
                ratio_str = np.array2string(ratio, formatter=formatter)+'%'
                Values_str = ' Values{:}, Ratios'.format(item_numbers)
                item_value = Values_str+ratio_str
        self.information_dict[item_name] = item_value

def _check_rainfall_rate_values(rain_source, times_in_1st_col=True):
    """ Check the rainfall rate values in rain source array
    rain_source: (numpy array) rainfall source array
          The 1st column is usually time series in seconds, from the 2nd column
          towards end columns are rainfall rate in m/s
    times_in_1st_col: indicate whether the first column is times
    Return:
        values_max: maximum rainfall rate in mm/h
        values_mean: average rainfall rate in mm/h
    """
    # get the pure rainfall rate values
    if times_in_1st_col:
        rain_values = rain_source[:, 1:]
        time_series = rain_source[:, 0]
    else:
        rain_values = rain_source
        time_series = np.arange(rain_values.shape[0])
    # convert the unit of rain rate values from m/s to mm/h
    rain_values_mmh = rain_values*3600*1000
    values_max = rain_values_mmh.max()
    values_mean = rain_values.mean(axis=1)
    rain_total_amount = np.trapz(y=values_mean, x=time_series) # in meter
    duration = np.ptp(time_series)
    rain_rate_mean = rain_total_amount*1000/(duration/3600) #mm/h
    if values_max > 5000 or rain_rate_mean > 1000:
        warnings.warn('Very large rainfall rates, better check your data!')
    return values_max, rain_rate_mean