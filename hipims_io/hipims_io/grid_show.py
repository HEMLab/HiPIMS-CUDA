#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
grid_show.py
To do:
    To visulize grid data, e.g. raster objec(s)
    static map functions:
    1. mapshow: general map of the grid values
    2. rankshow: show grid values in ranks
    3. hillshade: show a hillshade map of a grid
    4. vectorshow: show a vector map of two grids       
    animation functions:
    5. make_gif: create a gif file to show values of a series of grids
    6. make_mp4: create a video file to show values of a series of grids
Created on Tue Mar 10 15:37:28 2020

@author: Xiaodong Ming
"""
import os
import copy
import imageio
import numpy as np
import matplotlib.pyplot as plt
#plt.switch_backend('agg')
import matplotlib.colors as colors
import hipims_io.spatial_analysis as sp
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.colors import LightSource
from mpl_toolkits.axes_grid1 import make_axes_locatable
#%% draw inundation map with domain outline
def mapshow(raster_obj=None, array=None, header=None, ax=None,
            figname=None, figsize=None, dpi=300, 
            cax=True, cax_str=None, relocate=False, scale_ratio=1, **kwargs):
    """
    Display raster data without projection
    raster_obj: a Raster object
    array, header: to make Raster object if raster_obj is not given
    figname: the file name to export map, if figname is empty, then
        the figure will not be saved
    figsize: the size of map
    dpi: The resolution in dots per inch
    vmin and vmax define the data range that the colormap covers
    **kwargs: keywords argument of function imshow
    """
    if raster_obj is not None:
        array = raster_obj.array
        header = raster_obj.header
    # change NODATA_value to nan
    np.warnings.filterwarnings('ignore')
    array = array+0
    ind = array == header['NODATA_value']
    if ind.sum()>0:
        array = array.astype('float32')
        array[ind] = np.nan
    # adjust tick label and axis label
    map_extent = sp.header2extent(header)    
    if ax is None:
        fig, ax = plt.subplots(1, figsize=figsize)
    else:
        fig = ax.get_figure()
    img = plt.imshow(array, extent=map_extent, **kwargs)
    _adjust_axis_tick(ax, relocate, scale_ratio)
    # add colorbar
	# create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    if cax==True:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        cbar = plt.colorbar(img, cax=cax)
        plt.yticks(fontsize='small')
        if cax_str is not None:
            cbar.ax.set_xlabel(cax_str, horizontalalignment='left',
                               fontsize='small')
            cbar.ax.xaxis.set_label_coords(0, 1.06)
    ax.axes.grid(linestyle='-.', linewidth=0.2)
    ax.set_aspect('equal', 'box')
    # save figure
    if figname is not None:
        fig.savefig(figname, dpi=dpi)
    # return fig and axis handles
    return fig, ax

def rankshow(raster_obj=None, array=None, header=None, 
             breaks=[0.2, 0.3, 0.5, 1, 2], color='Blues',
             show_colorbar=True, show_colorlegend=False,
             figname=None, figsize=None, dpi=200,
             relocate=False, scale_ratio=1, **kwargs):
    """ Display water depth map in ranks defined by breaks
    color: color series of the ranks 
    """
    if raster_obj is not None:
        array = raster_obj.array
        header = raster_obj.header
    # change nan values to the min grid value
    array = array+0
    np.warnings.filterwarnings('ignore')
    ind = array == header['NODATA_value']
    array = array.astype('float32')
    array[ind] = np.nan
    min_value = np.nanmin(array)
    max_value = np.nanmax(array)
    array[np.isnan(array)] = 0
    # create color ranks
    breaks = copy.deepcopy(breaks)
    if breaks[0] > min_value:
        breaks.insert(0, min_value)
    if breaks[-1] < max_value:
        breaks.append(max_value)        
    norm = colors.BoundaryNorm(breaks, len(breaks))
    newcolors = cm.get_cmap(color, norm.N)
    newcolors = newcolors(np.linspace(0, 1, norm.N))
    newcolors[0, :] = np.array([255/256, 255/256, 255/256, 1]) #white for 1st
    newcmp = ListedColormap(newcolors)     
    map_extent = sp.header2extent(header)
    fig, ax = plt.subplots(figsize=figsize)
    chm_plot = ax.imshow(array, extent=map_extent, 
                         cmap=newcmp, norm=norm, alpha=0.7)
    _adjust_axis_tick(ax, relocate, scale_ratio)
    # create colorbar
    if show_colorbar is True:
        _set_colorbar(ax, chm_plot, norm)
    if show_colorlegend is True: # legend
        _set_color_legend(ax, norm, newcmp)
    if figname is not None:
        fig.savefig(figname, dpi=dpi)
    return fig, ax

def hillshade(raster_obj, figsize=None, azdeg=315, altdeg=45, vert_exag=1,
              alpha=1):
    """ Draw a hillshade map
    """
    array = raster_obj.array+0
    array[np.isnan(array)]=0
    ls = LightSource(azdeg=azdeg, altdeg=altdeg)
    cmap = plt.cm.gist_earth
    fig, ax = plt.subplots(figsize=figsize)
    rgb = ls.shade(array, cmap=cmap, 
                   blend_mode='overlay',vert_exag=vert_exag)
    ax.imshow(rgb, alpha=alpha)
    ax.set_axis_off()
    return fig, ax

def vectorshow(obj_x, obj_y, figname=None, figsize=None, dpi=300, **kwargs):
    """
    plot velocity map of U and V, whose values stored in two raster
    objects seperately
    """
    X, Y = obj_x.GetXYcoordinate()        
    U = obj_x.array
    V = obj_y.array
    if U.shape!=V.shape:
        raise TypeError('bad argument: the shapes must be the same')
    if 'figsize' in kwargs:
        figsize = kwargs['figsize']
    else:
        figsize = None
    fig, ax = plt.subplots(1, figsize=figsize)
    plt.quiver(X, Y, U, V)
    ax.set_aspect('equal', 'box')
    ax.tick_params(axis='y', labelrotation=90)
    if figname is not None:
        fig.savefig(figname, dpi=dpi)
    return fig, ax

def make_gif(output_file, obj_list=None, header=None, array_3d=None, 
                     time_str=None, breaks=None,
                     duration=0.5, **kwargs):
    """ Create animation of gridded data    
    mask_header: (dict) header file provide georeference of rainfall mask
    start_date: a datetime object to give the initial date and time of rain
    duration: duration for each frame (seconds)
    cellsize: sclar (meter) the size of rainfall grid cells
    """
    fig_names = _plot_temp_figs(obj_list, header, array_3d, breaks,
                                time_str, **kwargs)
    # create animation with the images
    images = []
    for fig_name in fig_names:
        images.append(imageio.imread(fig_name))
        os.remove(fig_name)
    # save animation and delete images
    if not output_file.endswith('.gif'):
        output_file = output_file+'.gif'
    imageio.mimsave(output_file, images, duration=duration)

def make_mp4(output_file, obj_list=None, header=None, array_3d=None, 
               time_str=None, breaks=None, fig_names=None,
               fps=10, **kwargs):
    """ Create a video file based on a series of grids
    obj_list: a list of Raster objects
    header: a header dict providing georeference the grid [not necessary if 
                                                           obj_list was given]
    array_3d: a 3D numpy array storing grid values for each timestep 
                (in 1st dimension), [not necessary if obj_list was given]
    time_str: a list of string to show time information for each frame
    """
    if fig_names is None:
        fig_names = _plot_temp_figs(obj_list, header, array_3d, breaks,
                                    time_str, **kwargs)
    if not output_file.endswith('.mp4'):
        output_file = output_file+'.mp4'
    print(output_file)
    writer = imageio.get_writer(output_file, 'MP4', fps=fps)
    for fig_name in fig_names:
        writer.append_data(imageio.imread(fig_name))
        os.remove(fig_name)
    writer.close()

def plot_shape_file(shp_file, color='r', linewidth=0.5, figsize=None, ax=None):
    """plot a shape file to a map axis
    """
    import shapefile
    sf = shapefile.Reader(shp_file)
    if ax is None:
        fig, ax = plt.subplots(1, figsize=figsize)
        xbound = None
        ybound = None
    else:
        fig = ax.get_figure()
        xbound = ax.get_xbound()
        ybound = ax.get_ybound()
    # draw shape file on the rainfall map
    for shape in sf.shapeRecords():
        for i in range(len(shape.shape.parts)):
            i_start = shape.shape.parts[i]
            if i==len(shape.shape.parts)-1:
                i_end = len(shape.shape.points)
            else:
                i_end = shape.shape.parts[i+1]
            x = [i[0] for i in shape.shape.points[i_start:i_end]]
            y = [i[1] for i in shape.shape.points[i_start:i_end]]
            ax.plot(x, y, color=color, linewidth=linewidth)
    ax.set_xbound(xbound)
    ax.set_ybound(ybound)
    return fig, ax

def _plot_temp_figs(obj_list=None, header=None, array_3d=None,
                    breaks=None, time_str=None, **kwargs):
    """plot a series of temp pictures and save to make animation
    plot_fun: the function to plot
    """
    """ Create a video file based on a series of grids
    """
    if obj_list is not None:
        header = obj_list[0].header
        array_3d = []
        for grid_obj in obj_list:
            array_3d.append(grid_obj.array)
        array_3d = np.array(array_3d)
    if breaks is None: #plot continous map
        plot_fun = mapshow
    else:
        plot_fun = rankshow
        kwargs['breaks']=breaks
    fig_names = []
    for i in np.arange(array_3d.shape[0]):
        fig_name = 'temp'+str(i)+'.png'
        fig, ax = plot_fun(array=array_3d[i],header=header, **kwargs)
        if type(time_str) is list:
            ax.set_title(time_str[i])
        fig.savefig(fig_name)
        plt.close(fig)
        fig_names.append(fig_name)
    return fig_names
    
#%%
def _set_colorbar(ax, img, norm):
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

def _adjust_axis_tick(ax, relocate=True, scale_ratio=1):
    """
    Adjust the axis tick to a new staring point and/or new unit 
    Example:
        if scale_ratio = 1000, and the original extent unit is meter,
        then the unit is converted to km, and the extent is divided by 1000
    """
    xticks = ax.get_xticks()
    x_space = xticks[1]-xticks[0]
    yticks = ax.get_yticks()
    y_space = yticks[1]-yticks[0]
    if relocate:
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        xticks = np.arange(xlim[0], xlim[1], x_space)
        xticks_label = xticks-xlim[0]
        yticks = np.arange(ylim[0], ylim[1], y_space)
        yticks_label = yticks-ylim[0]
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
    else:
        xticks_label = xticks
        yticks_label = yticks
    if scale_ratio == 1000:
        label_tag = 'km'
        xticks_label = (xticks_label/1000).astype('int64')
        yticks_label = (yticks_label/1000).astype('int64')
    else:
        label_tag = 'meter'
        xticks_label = (xticks_label).astype('int64')
        yticks_label = (yticks_label).astype('int64')
    ax.set_xticklabels(xticks_label)
    ax.set_yticklabels(yticks_label)
    ax.set_xlabel(label_tag+' towards east')
    ax.set_ylabel(label_tag+' towards north')
    return None

def main():
    print('Package to show grid data')

if __name__=='__main__':
    main()