#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
spatial_analysis
functions to analyse data in raster and/or feature datasets
to replace ArcGridDataProcessing
-------------------------------------------------------------------------------
@author: Xiaodong Ming
Created on Wed Nov  6 14:33:36 2019
-------------------------------------------------------------------------------
Assumptions:
- map unit is meter
- its cellsize is the same in both x and y direction
- its reference position is on the lower left corner of the southwest cell
To do:
- read and write arc
"""
__author__ = "Xiaodong Ming"
import os
import gzip
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.colors import LightSource
from mpl_toolkits.axes_grid1 import make_axes_locatable
#%%
def arc_header_read(file_name, header_rows=6):
    """ read the header of a asc file as a dictionary
    file_name: (string) file name
    header_rows: (int) number of header rows
    Return:
        header: a dictionary with keys:
            ncols: (int) number of columns
            nrows: (int) number of rows
            xllcorner: (int/float) x-coordinate of the lower left corner of
                the lower left cell of the grid
            yllcorner: (int/float) y-coordinate of the lower left corner of
                the bottom left cell of the grid
            cellsize: (int/float) the length of one square cell
            NODATA_value: (int/float)|-9999 the value representing nodata cell
    """
    check_file_existence(file_name)
    # read header
    header = {} # store header information including ncols, nrows,...
    row_ite = 1
    if file_name.endswith('.gz'):
        with gzip.open(file_name, 'rt') as file_h:     
            for line in file_h:
                if row_ite <= header_rows:
                    line = line.split(" ", 1)
                    header[line[0]] = float(line[1])
                else:
                    break
                row_ite = row_ite+1
    else:
        # read header
        with open(file_name, 'rt') as file_h:
            for line in file_h:
                if row_ite <= header_rows:
                    line = line.split(" ", 1)
                    header[line[0]] = float(line[1])
                else:
                    break
                row_ite = row_ite+1
    header['ncols'] = int(header['ncols'])
    header['nrows'] = int(header['nrows'])
    return header

def arcgridread(file_name, header_rows=6, return_nan=True):
    """ Read ArcGrid format raster file
    file_name: (str) the file name to read data
    header_rows: (int) the number of head rows of the asc file
    Return:
        array: (int/float numpy array)
        header: (dict) to provide reference information of the grid
        extent: (tuple) outline extent of the grid (left, right, bottom, top)
    Note: this function can also read compressed gz files
    """
    check_file_existence(file_name)
    # read header
    header = arc_header_read(file_name, header_rows)
    # read value array
    array = np.loadtxt(file_name, skiprows=header_rows, dtype='float64')
    if return_nan:
        array[array == header['NODATA_value']] = np.nan
    prj_file = file_name[:-4]+'.prj'
    if os.path.isfile(prj_file):
        with open(prj_file, 'r') as file:
            projection = file.read()
    else:
        projection = None
    return array, header, projection

def arcgridwrite(file_name, array, header, compression=False):
    """ write gird data into a ascii file
    file_name: (str) the file name to write grid data. A compressed file will
        automatically add a suffix '.gz'
    array: (int/float numpy array)
    header: (dict) to provide reference information of the grid
    compression: (logic) to inidcate whether compress the ascii file
    Example:
        gird = np.zeros((5,10))
        grid[0,:] = -9999
        grid[-1,:] = -9999
        header = {'ncols':10, 'nrows':5, 'xllcorner':0, 'yllcorner':0,
                  'cellsize':2, 'NODATA_value':-9999}
        file_name = 'example_file.asc'
        arcgridwrite(file_name, array, header, compression=False)
        arcgridwrite(file_name, array, header, compression=True)
    """
    array = array+0
    if not isinstance(header, dict):
        raise TypeError('bad argument: head')
    if file_name.endswith('.gz'):
        compression = True
    # create a file (gz or asc)
    if compression:
        if not file_name.endswith('.gz'):
            file_name = file_name+'.gz'
        file_h = gzip.open(file_name, 'wb')
    else:
        file_h = open(file_name, 'wb')
    file_h.write(b"ncols    %d\n" % header['ncols'])
    file_h.write(b"nrows    %d\n" % header['nrows'])
    file_h.write(b"xllcorner    %g\n" % header['xllcorner'])
    file_h.write(b"yllcorner    %g\n" % header['yllcorner'])
    file_h.write(b"cellsize    %g\n" % header['cellsize'])
    file_h.write(b"NODATA_value    %g\n" % header['NODATA_value'])
    array[np.isnan(array)] = header['NODATA_value']
    np.savetxt(file_h, array, fmt='%g', delimiter=' ')
    file_h.close()
    print(file_name + ' created')

def tif_read(file_name):
    """
    read tif file and return array, header, projection
    only read the first band
    """
    from osgeo import gdal
    ds = gdal.Open(file_name)        
    ncols = ds.RasterXSize
    nrows = ds.RasterYSize
    geo_transform = ds.GetGeoTransform()
    x_min = geo_transform[0]
    cellsize = geo_transform[1]
    y_max = geo_transform[3]
    xllcorner = x_min
    yllcorner = y_max-nrows*cellsize
    rasterBand = ds.GetRasterBand(1)
    NODATA_value = rasterBand.GetNoDataValue()
    array = rasterBand.ReadAsArray()
    header = {'ncols':ncols, 'nrows':nrows,
              'xllcorner':xllcorner, 'yllcorner':yllcorner,
              'cellsize':cellsize, 'NODATA_value':NODATA_value}     
    if not np.isscalar(header['NODATA_value']):
        header['NODATA_value'] = -9999
    array[array == header['NODATA_value']] = float('nan')
    projection = ds.GetProjection()
    rasterBand = None
    ds = None
    return array, header, projection

def byte_file_read(file_name):
    """ Read file from a bytes object
    """
    # read header
    header = {} # store header information including ncols, nrows, ...
    num_header_rows = 6
    for _ in range(num_header_rows):
        line = file_name.readline()
        line = line.strip().decode("utf-8").split(" ", 1)
        header[line[0]] = float(line[1])
        # read value array
    array  = np.loadtxt(file_name, skiprows=num_header_rows, 
                        dtype='float64')
    array[array == header['NODATA_value']] = float('nan')
    header['ncols'] = int(header['ncols'])
    header['nrows'] = int(header['nrows'])
    return array, header

#%% Combine raster files
def combine_raster(asc_files, num_header_rows=6):
    """Combine a list of asc files to a DEM Raster
    asc_files: a list of asc file names
    all raster files have the same cellsize
    """
    # default values for the combined Raster file
    xllcorner_all = []
    yllcorner_all = []
    extent_all =[]
    # read header
    for file in asc_files:
        header0 = arc_header_read(file, num_header_rows)
        extent0 = header2extent(header0)
        xllcorner_all.append(header0['xllcorner'])
        yllcorner_all.append(header0['yllcorner'])
        extent_all.append(extent0)
    cellsize = header0['cellsize']
    if 'NODATA_value' in header0.keys():
        NODATA_value = header0['NODATA_value']
    else:
        NODATA_value = -9999
    xllcorner_all = np.array(xllcorner_all)
    xllcorner = xllcorner_all.min()
    yllcorner_all = np.array(yllcorner_all)
    yllcorner = yllcorner_all.min()
    extent_all = np.array(extent_all)
    x_min = np.min(extent_all[:,0])
    x_max = np.max(extent_all[:,1])
    y_min = np.min(extent_all[:,2])
    y_max = np.max(extent_all[:,3])
#    extent = (x_min, x_max, y_min, y_max)
#    print(extent)
    nrows = int((y_max-y_min)/cellsize)
    ncols = int((x_max-x_min)/cellsize)
    header = header0.copy()
    header['xllcorner'] = xllcorner
    header['yllcorner'] = yllcorner
    header['ncols'] = ncols
    header['nrows'] = nrows
    header['NODATA_value'] = NODATA_value
    array = np.zeros((nrows ,ncols))+NODATA_value
    print(array.shape)
    for file in asc_files:
        array0, header0, _ = arcgridread(file, num_header_rows)
        extent0 = header2extent(header0)
        x0 = extent0[0]+header0['cellsize']/2
        y0 = extent0[3]-header0['cellsize']/2
        row0, col0 = map2sub(x0, y0, header)
        array[row0:row0+header0['nrows'],
              col0:col0+header0['ncols']] = array0
    array[array == header['NODATA_value']] = float('nan')
    extent = header2extent(header)
    return array, header, extent


#%% ----------------------------Visulization-----------------------------------
def map_show(array, header, figname=None, figsize=None, dpi=300,
             vmin=None, vmax=None,
             cax=True, relocate=False, scale_ratio=1):
    """
    Display raster data
    figname: the file name to export map, if figname is empty, then
        the figure will not be saved
    figsize: the size of map
    dpi: The resolution in dots per inch
    vmin and vmax define the data range that the colormap covers
    """
    np.warnings.filterwarnings('ignore')
    array = array+0
    fig, ax = plt.subplots(1, figsize=figsize)
    # draw grid data
    array[array==header['NODATA_value']]=np.nan
    # adjust tick label and axis label
    map_extent = header2extent(header)
    map_extent = _adjust_map_extent(map_extent, relocate, scale_ratio)
    img=plt.imshow(array, extent=map_extent, vmin=vmin, vmax=vmax)
    # colorbar
	# create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    if cax==True:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(img, cax=cax)
    ax.axes.grid(linestyle='-.', linewidth=0.2)
    # save figure
    if figname is not None:
        fig.savefig(figname, dpi=dpi)
    return fig, ax

def rank_show(array, header, figname=None, figsize=None, dpi=300,
            breaks=[0.2, 0.3, 0.5, 1, 2], # default for water depth
            show_colorbar=True, show_colorlegend=False,
            relocate=False, scale_ratio=1):
    """ 
    Categorize array data as ranks according to the breaks and display a ranked
        map
    """
    np.warnings.filterwarnings('ignore')
    array = array+0
    if 'NODATA_value' in header.keys():
        array[array == header['NODATA_value']] = np.nan
    if breaks[0] > np.nanmin(array):
        breaks.insert(0, np.nanmin(array))
    if breaks[-1] < np.nanmax(array):
        breaks.append(np.nanmax(array))        
    norm = colors.BoundaryNorm(breaks, len(breaks))
    blues = cm.get_cmap('blues', norm.N)
    newcolors = blues(np.linspace(0, 1, norm.N))
    white = np.array([255/256, 255/256, 255/256, 1])
    newcolors[0, :] = white
    newcmp = ListedColormap(newcolors)
    map_extent = header2extent(header)
    map_extent = _adjust_map_extent(map_extent, relocate, scale_ratio)
    fig, ax = plt.subplots(figsize=figsize)
    chm_plot = ax.imshow(array, extent=map_extent, 
                         cmap=newcmp, norm=norm, alpha=0.7)
    # create colorbar
    if show_colorbar is True:
        _set_colorbar(ax, chm_plot, norm)
    if show_colorlegend is True: # legend
        _set_color_legend(ax, norm, newcmp)
    plt.show()
    # save figure
    if figname is not None:
        fig.savefig(figname, dpi=dpi)
    return fig, ax

def hillshade_show(array, header, figsize=None,
                   azdeg=315, altdeg=45, vert_exag=1):
    """ Draw a hillshade map
    """
    array = array+0
    array[np.isnan(array)] = 0
    array[array == header['NODATA_value']] = 0
    map_extent = header2extent(header)
    ls = LightSource(azdeg=azdeg, altdeg=altdeg)
    cmap = plt.cm.gist_earth
    fig, ax = plt.subplots(figsize=figsize)
    rgb = ls.shade(array, cmap=cmap, 
                   blend_mode='overlay',vert_exag=vert_exag)
    ax.imshow(rgb, extent=map_extent)
#    ax.set_axis_off()
    plt.show()
    return fig, ax

#%% ------------------------Supporting functions-------------------------------
def check_file_existence(file_name):
    """ check whether a file exists
    """
    try:
        file_h = open(file_name, 'r')
        file_h.close()
    except FileNotFoundError:
        raise

def header2extent(header):
    """ convert a header dict to a 4-element tuple (left, right, bottom, top)
    all four elements shows the coordinates at the edge of a cell, not center
    """
    left = header['xllcorner']
    right = header['xllcorner']+header['ncols']*header['cellsize']
    bottom = header['yllcorner']
    top = header['yllcorner']+header['nrows']*header['cellsize']
    extent = (left, right, bottom, top)
    return extent

def shape_extent_to_header(shape, extent, nan_value=-9999):
    """ Create a header dict with shape and extent of an array
    """
    ncols = shape[1]
    nrows = shape[0]
    xllcorner = extent[0]
    yllcorner = extent[2]
    cellsize_x = (extent[1]-extent[0])/ncols
    cellsize_y = (extent[3]-extent[2])/nrows
    if cellsize_x != cellsize_y:
        raise ValueError('extent produces different cellsize in x and y')
    cellsize = cellsize_x
    header = {'ncols':ncols, 'nrows':nrows,
              'xllcorner':xllcorner, 'yllcorner':yllcorner,
              'cellsize':cellsize, 'NODATA_value':nan_value}
    return header

def map2sub(X, Y, header):
    """ convert map coordinates to subscripts of an array
    array is defined by a geo-reference header
    X, Y: a scalar or numpy array of coordinate values
    Return: rows, cols in the array
    """
    # X and Y coordinate of the centre of the first cell in the array
    x0 = header['xllcorner']+0.5*header['cellsize']
    y0 = header['yllcorner']+(header['nrows']-0.5)*header['cellsize']
    rows = (y0-Y)/header['cellsize'] # row and col number starts from 0
    cols = (X-x0)/header['cellsize']
    if isinstance(rows, np.ndarray):
        rows = rows.astype('int64')
        cols = cols.astype('int64') #.astype('int64')
    else:
        rows = int(rows)
        cols = int(cols)
    return rows, cols

def sub2map(rows, cols, header):
    """
    convert subscripts of a matrix to map coordinates 
    rows, cols: subscripts of the data matrix, starting from 0
    return
        X, Y: coordinates in map units
    """
    #x and y coordinate of the centre of the first cell in the matrix
    if not isinstance(rows, np.ndarray):
        rows = np.array(rows)
        cols = np.array(cols)        
    extent = header2extent(header)
    left = extent[0] #(left, right, bottom, top)
    top = extent[3]
    X = left + (cols+0.5)*header['cellsize']
    Y = top  - (rows+0.5)*header['cellsize']  
    return X, Y

#% Extent compare between two Raster objects
def compare_extent(extent0, extent1):
    """Compare and show the difference between two Raster extents
    extent0, extent1: objects or extent dicts to be compared
    displaye: whether to show the extent in figures
    Return:
        0 extent0>=extent1
        1 extent0<extent1
        2 extent0 and extent1 have intersections
    """
    logic_left = extent0[0]<=extent1[0]
    logic_right = extent0[1]>=extent1[1]
    logic_bottom = extent0[2]<=extent1[2]
    logic_top = extent0[3]>=extent1[3]
    logic_all = logic_left+logic_right+logic_bottom+logic_top
    if logic_all == 4:
        output = 0
    elif logic_all == 0:
        output = 1
    else:
        output = 2
        print(extent0)
        print(extent1)
    return output
   
def extent2shape_points(extent):
    """Convert extent to a two-col numpy array of shape points
    """
    #extent = (left, right, bottom, top)
    shape_points = np.array([[extent[0], extent[2]], 
                             [extent[1], extent[2]], 
                             [extent[1], extent[3]], 
                             [extent[0], extent[3]]])
    return shape_points

def _adjust_map_extent(extent, relocate=True, scale_ratio=1):
    """
    Adjust the extent (left, right, bottom, top) to a new staring point 
        and new unit. extent values will be divided by the scale_ratio
    Example:
        if scale_ratio = 1000, and the original extent unit is meter,
        then the unit is converted to km, and the extent is divided by 1000
    """
    if relocate:
        left = 0 
        right = (extent[1]-extent[0])/scale_ratio
        bottom = 0
        top = (extent[3]-extent[2])/scale_ratio
    else:
        left = extent[0]/scale_ratio
        right = extent[1]/scale_ratio
        bottom = extent[2]/scale_ratio
        top = extent[3]/scale_ratio
    return (left, right, bottom, top)

def _set_colorbar(ax,img,norm):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(img, cax=cax)
    y_tick_values = cax.get_yticks()
    boundary_means = [np.mean((y_tick_values[ii],y_tick_values[ii-1])) 
                        for ii in range(1, len(y_tick_values))]
    category_names = [(str(norm.boundaries[ii-1])+'~'+
                       str(norm.boundaries[ii]))
                      for ii in range(1, len(norm.boundaries))]
    category_names[0] = '<='+str(norm.boundaries[1])
    category_names[-1] = '>'+str(norm.boundaries[-2])
    cax.yaxis.set_ticks(boundary_means)
    cax.yaxis.set_ticklabels(category_names,rotation=0)
    return cax

def _set_color_legend(ax, norm, cmp,
                      loc='lower right', bbox_to_anchor=(1,0),
                      facecolor=None):
    category_names = [(str(norm.boundaries[ii-1])+'~'+
                       str(norm.boundaries[ii]))
                      for ii in range(1, len(norm.boundaries))]
    category_names[0] = '<='+str(norm.boundaries[1])
    category_names[-1] = '>'+str(norm.boundaries[-2])
    ii = 0
    legend_labels = {}
    for category_name in category_names:
        legend_labels[category_name] = cmp.colors[ii,]
        ii = ii+1
    patches = [Patch(color=color, label=label)
               for label, color in legend_labels.items()]
    ax.legend(handles=patches, loc=loc,
              bbox_to_anchor=bbox_to_anchor,
              facecolor=facecolor)
    return ax

def main():
    print('Fucntions to process asc data')

if __name__=='__main__':
    main()