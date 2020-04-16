#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
raster
To do:
    Read, write and analyze gridded Raster data
Created on Tue Mar 31 16:20:13 2020

@author: Xiaodong Ming
"""
import copy
import math
import numpy as np
#from osgeo import gdal, ogr, osr
from scipy import interpolate
import hipims_io.spatial_analysis as sp
import hipims_io.grid_show as gs
#%% *******************************To deal with raster data********************
#   ***************************************************************************    
class Raster(object):
    """    
    To deal with raster data with a ESRI ASCII or GTiff format
    Properties:
        source_file: file name to read grid data
        output_file: file name to write a raster object
        array: a numpy array storing grid cell values
        header: a dict storing reference information of the grid
        extent: a tuple storing outline limits of the raster (left, right, 
        bottom, top)
        extent_dict: a dictionary storing outline limits of the raster
        projection: (string) the Well-Known_Text (wkt) projection information
    
    Methods(public):
        Write_asc: write grid data into an asc file with or without 
            compression(.gz)
        to_osgeo_raster: convert this object to an osgeo raster object
        rect_clip: clip raster according to a rectangle extent
        clip: clip raster according to a polygon
        rasterize: rasterize a shapefile on the Raster object and return a 
            bool array with 'Ture' in and on the polygon/polyline
        resample: resample the raster to a new cellsize
        GetXYcoordinate: Get X and Y coordinates of all raster cells
        mapshow: draw a map of the raster dataset
        VelocityShow: draw velocity vectors as arrows with values on two Raster
            datasets (u, v)
            
    Methods(private):
        __header2extent: convert header to extent
        __read_asc: read an asc file ends with .asc or .gz
            with a reference header
        __read_tif: read tiff file
        
    """
#%%======================== initialization function ===========================   
    def __init__(self, source_file=None, array=None, header=None, 
                 epsg=None, projection=None, num_header_rows=6):
        """
        source_file: name of a asc/tif file if a file read is needed
        array: values in each raster cell [a numpy array]
        header: georeference of the raster [a dictionary containing 6 keys]:
            nrows, nclos [int]
            cellsize, xllcorner, yllcorner
            NODATA_value
        epsg: epsg code [int]
        projection: WktProjection [string]
        """
        if epsg is not None:
            projection = self.__set_wkt_projection(epsg)
        if type(source_file) is str:
            if source_file.endswith('.tif'):
                array, header, projection = sp.tif_read(source_file) # only read the first band
            else:
                array, header, projection = sp.arcgridread(source_file,
                                                    num_header_rows)
            self.source_file = source_file
        elif type(source_file) is bytes:  # try a binary file-like object
            array, header = sp.byte_file_read(source_file)
        extent = sp.header2extent(header)
        self.source_file = source_file
        self.projection = projection
        self.array = array
        self.header = header
        self.extent = extent
        self.extent_dict = {'left':extent[0], 'right':extent[1],
                            'bottom':extent[2], 'top':extent[3]}
            
#%%============================= Spatial analyst ==============================   
    def rect_clip(self, clip_extent):
        """
        clip_extent: left, right, bottom, top
        clip raster according to a rectangle extent
        return:
           a new raster object
        """
        X = clip_extent[0:2]
        Y = clip_extent[2:4]
        rows, cols = sp.map2sub(X, Y, self.header)
        x_centre, y_centre = sp.sub2map(rows, cols, self.header)
        xllcorner = min(x_centre)-0.5*self.header['cellsize']
        yllcorner = min(y_centre)-0.5*self.header['cellsize']
        header_new = copy.deepcopy(self.header)
        array_new = self.array[min(rows):max(rows), min(cols):max(cols)]
        header_new['nrows'] = array_new.shape[0]
        header_new['ncols'] = array_new.shape[1]
        header_new['xllcorner'] = xllcorner
        header_new['yllcorner'] = yllcorner
        new_obj = Raster(array=array_new, header=header_new,
                         projection=self.projection)
        return new_obj
    
    def clip(self, mask=None):
        """
        clip raster according to a mask
        mask: 
            1. string name of a shapefile
            2. numpy vector giving X and Y coords of the mask points
        
        return:
            a new raster object
        """
        if isinstance(mask, str):
            shpName =  mask
        # Open shapefile datasets  
        from osgeo import ogr
        shpDriver = ogr.GetDriverByName('ESRI Shapefile')
        shpDataset = shpDriver.Open(shpName, 0) # 0=Read-only, 1=Read-Write
        layer = shpDataset.GetLayer()
        shpExtent = np.array(layer.GetExtent()) #(minX, maxY, maxX, minY)           
        # 1. rectangle clip raster
        new_obj = self.rect_clip(shpExtent)
        new_raster = copy.deepcopy(new_obj)                
        indexArray = new_raster.rasterize(shpDataset)
        arrayClip = new_raster.array
        arrayClip[indexArray==0]=new_raster.header['NODATA_value']
        new_raster.array = arrayClip        
        shpDataset.Destroy()
        return new_raster
    
    def rasterize(self, shpDSName, rasterDS=None):
        """
        rasterize the shapefile to the raster object and return a bool array
            with Ture value in and on the polygon/polyline
        shpDSName: string for shapefilename, dataset for ogr('ESRI Shapefile')
            object
        return numpy array
        """
        from osgeo import gdal, ogr
        if isinstance(shpDSName, str):
            shpDataset = ogr.Open(shpDSName)
        else:
            shpDataset = shpDSName
        layer = shpDataset.GetLayer()
        if rasterDS is None:
            obj_raster = copy.deepcopy(self)
            obj_raster.array = np.zeros(obj_raster.array.shape)
            target_ds = obj_raster.to_osgeo_raster()
        else:
            target_ds = rasterDS
        gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[-1])
        rasterized_array = target_ds.ReadAsArray()
        indexArray = np.full(rasterized_array.shape, False)
        indexArray[rasterized_array==-1] = True
        target_ds=None
        return indexArray
    
    def resample(self, cellsize_n, method='bilinear'):
        """
        resample the raster to a new cellsize
        cellsize_n: cellsize of the new raster
        method: Resampling method to use. Available methods are:
            near: nearest neighbour resampling (default, fastest algorithm, 
                                                worst interpolation quality).        
            bilinear: bilinear resampling.        
            cubic: cubic resampling.        
            cubicspline: cubic spline resampling.        
            lanczos: Lanczos windowed sinc resampling.        
            average: average resampling, computes the average of all 
                    non-NODATA contributing pixels.        
            mode: mode resampling, selects the value which appears most often 
                    of all the sampled points.        
            max: maximum resampling, selects the maximum value from all 
                    non-NODATA contributing pixels.        
            min: minimum resampling, selects the minimum value from all 
                    non-NODATA contributing pixels.        
            med: median resampling, selects the median value of all 
                    non-NODATA contributing pixels.        
            q1: first quartile resampling, selects the first quartile 
                value of all non-NODATA contributing pixels.        
            q3: third quartile resampling, selects the third quartile 
                value of all non-NODATA contributing pixels
        """
        from osgeo import gdal
        cellSize = self.header['cellsize']
        ras_x_size = self.header['ncols']
        newras_x_size = int(ras_x_size*cellSize/cellsize_n)
        rasterYSize = self.header['nrows']
        newRasterYSize = int(rasterYSize*cellSize/cellsize_n)
        
        g = self.to_osgeo_raster() # get original gdal dataset
        total_obs = g.RasterCount
        drv = gdal.GetDriverByName( "MEM" )
        dst_ds = drv.Create('', g.RasterXSize, g.RasterYSize, 1,
                            eType=gdal.GDT_Float32)
        dst_ds.SetGeoTransform( g.GetGeoTransform())
        dst_ds.SetProjection ( g.GetProjectionRef() )
        hires_data = self.array
        dst_ds.GetRasterBand(1).WriteArray ( hires_data )
        
        geo_trans_v = g.GetGeoTransform()
        drv = gdal.GetDriverByName( "MEM" )
        resampled_ds = drv.Create('', newras_x_size, newRasterYSize, 1, 
                                  eType=gdal.GDT_Float32)

        geo_trans_v_new = (geo_trans_v[0], cellsize_n, geo_trans_v[2],
                           geo_trans_v[3], geo_trans_v[3], -cellsize_n)
        resampled_ds.SetGeoTransform(geo_trans_v_new )
        resampled_ds.SetProjection (g.GetProjectionRef() )
        resampled_ds.SetMetadata ({"TotalNObs":"%d" % total_obs})

        gdal.RegenerateOverviews(dst_ds.GetRasterBand(1),
                                 [resampled_ds.GetRasterBand(1)], method)
    
        resampled_ds.GetRasterBand(1).SetNoDataValue(self.header['NODATA_value'])
        
        new_obj = self.__osgeo2raster(resampled_ds)
        resampled_ds = None

        return new_obj
    
    def point_interpolate(self, points, values, method='nearest'):
        """ Interpolate values of 2D points to all cells on the Raster object
        2D interpolate
        points: ndarray of floats, shape (n, 2)
            Data point coordinates. Can either be an array of shape (n, 2), 
            or a tuple of ndim arrays.
        values: ndarray of float or complex, shape (n, )
            Data values.
        method: {‘linear’, ‘nearest’, ‘cubic’}, optional
            Method of interpolation.
        """
        grid_x, grid_y = self.GetXYcoordinate()
        array_interp = interpolate.griddata(points, values, (grid_x, grid_y),
                                            method=method)
        new_obj = copy.deepcopy(self)
        new_obj.array = array_interp
        new_obj.source_file = 'mask_'+new_obj.source_file
        return new_obj
    
    def grid_interpolate(self, value_grid, method='nearest'):
        """ Interpolate values of a grid to all cells on the Raster object
        2D interpolate
        value_grid: a grid file string or Raster object 
        method: {‘linear’, ‘nearest’, ‘cubic’}, optional
            Method of interpolation.
        Return: 
            a numpy array with the same size of the self object
        """
        if type(value_grid) is str:
            value_grid = Raster(value_grid)
        points_x, points_y = value_grid.GetXYcoordinate()
        points = np.c_[points_x.flatten(), points_y.flatten()]
        values = value_grid.array.flatten()
        ind_nan = ~np.isnan(values)
        grid_x, grid_y = self.GetXYcoordinate()
        array_interp = interpolate.griddata(points[ind_nan, :], values[ind_nan],
                                            (grid_x, grid_y), method=method)
        return array_interp
    
    def grid_resample(self, newsize):
        """
        resample a grid to a new grid resolution via nearest interpolation
        """
        if isinstance(newsize, dict):
            header = newsize.copy()
        else:
            oldsize = self.header['cellsize']
            header = copy.deepcopy(self.header)
            header['cellsize'] = newsize
            ncols = math.floor(oldsize*self.header['ncols']/newsize)
            nrows = math.floor(oldsize*self.header['nrows']/newsize)
            header['ncols'] = ncols
            header['nrows'] = nrows
        #centre of the first cell in array
        x11 = header['xllcorner']+0.5*header['cellsize']
        y11 = header['yllcorner']+(header['nrows']-0.5)*header['cellsize']
        x_all = np.linspace(x11, x11+(header['ncols']-1)*header['cellsize'],
                            header['ncols'])
        y_all = np.linspace(y11, y11-(header['nrows']-1)*header['cellsize'],
                            header['nrows'])
        row_all, col_all = sp.map2sub(x_all, y_all, self.header)
        rows, cols = np.meshgrid(row_all, col_all) # nrows*ncols array
        array = self.array[rows, cols]
        array = array.transpose()
        array = array.astype(self.array.dtype)
        new_obj = Raster(array=array, header=header)
        return new_obj
    
    def assign_to(self, new_header):
        """ Assign_to the object to a new grid defined by new_header 
        If cellsize are not equal, the origin Raster will be firstly 
        resampled to the target grid.
        obj_origin, obj_target: Raster objects
        """
        obj_origin = copy.deepcopy(self)
        if obj_origin.header['cellsize'] != new_header['cellsize']:
            obj_origin = obj_origin.GridResample(new_header['cellsize'])
        grid_x, grid_y = obj_origin.GetXYcoordinate()
        rows, cols = sp.map2sub(grid_x, grid_y, new_header)
        ind_r = np.logical_and(rows >= 0, rows <= new_header['nrows']-1)
        ind_c = np.logical_and(cols >= 0, cols <= new_header['ncols']-1)
        ind = np.logical_and(ind_r, ind_c)
#        ind = np.logical_and(ind, ~np.isnan(obj_origin.array))
        array = obj_origin.array[ind]
        array = np.reshape(array, (new_header['nrows'], new_header['ncols']))
#        array[rows[ind], cols[ind]] = obj_origin.array[ind]
        obj_output = Raster(array=array, header=new_header)
        return obj_output

    def to_points(self):
        """ Get X and Y coordinates of all raster cells
        return xv, yv numpy array with the same size of the raster object
        """
        ny, nx = self.array.shape
        cellsize = self.header['cellsize']
        # coordinate of the centre on the top-left pixel
        x00centre = self.extent_dict['left'] + cellsize/2
        y00centre = self.extent_dict['top'] - cellsize/2
        x = np.arange(x00centre, x00centre+cellsize*nx, cellsize)
        y = np.arange(y00centre, y00centre-cellsize*ny, -cellsize)
        xv, yv = np.meshgrid(x, y)
        return xv, yv
    
    def write_asc(self, output_file, EPSG=None, compression=False):
        
        """
        write raster as asc format file 
        output_file: output file name
        EPSG: epsg code, if it is given, a .prj file will be written
        compression: logic, whether compress write the asc file as gz
        """
        sp.arcgridwrite(output_file, self.array, self.header, compression)
        if EPSG is not None:
            self.__set_wkt_projection(EPSG)
        # if projection is defined, write .prj file for asc file
        if output_file.endswith('.asc'):
            if self.projection is not None:
                prj_file=output_file[0:-4]+'.prj'
                wkt = self.projection
                with open(prj_file, "w") as prj:        
                    prj.write(wkt)
        return None
    
    def to_osgeo_raster(self, filename=None, fileformat = 'GTiff',
                        destEPSG=27700):        
        """
        convert this object to an osgeo raster object, write a tif file if 
            necessary
        filename: the output file name
        fileformat: GTiff or AAIGrid
        destEPSG: the EPSG projection code default: British National Grid
        return:
            an osgeo raster dataset
            or a tif filename if it is written
        """
        from osgeo import gdal, osr
        if filename is None:
            dst_filename = ''
            driver_name = 'MEM'
        else:
            dst_filename = filename
            driver_name = fileformat
        if not dst_filename.endswith('.tif'):
            dst_filename = dst_filename+'.tif'
    
        # You need to get those values like you did.
        PIXEL_SIZE = self.header['cellsize']  # size of the pixel...        
        x_min = self.extent[0] # left  
        y_max = self.extent[3] # top
        dest_crs = osr.SpatialReference()
        dest_crs.ImportFromEPSG(destEPSG)
        # create dataset with driver
        driver = gdal.GetDriverByName(driver_name)
        ncols = int(self.header['ncols'])
        nrows = int(self.header['nrows'])
        dataset = driver.Create(dst_filename, 
            xsize=ncols, 
            ysize=nrows, 
            bands=1, 
            eType=gdal.GDT_Float32)
    
        dataset.SetGeoTransform((
            x_min,       # 0
            PIXEL_SIZE,  # 1
            0,           # 2
            y_max,       # 3
            0,           # 4
            -PIXEL_SIZE))  
    
        dataset.SetProjection(dest_crs.ExportToWkt())
        array = self.array
        dataset.GetRasterBand(1).WriteArray(array)
        dataset.GetRasterBand(1).SetNoDataValue(self.header['NODATA_value'])
        if filename is not None:
            dataset.FlushCache()  # Write to disk.
            dataset = None
            return dst_filename
        else:
            return dataset
#%%=============================Visualization==================================
    def mapshow(self, **kwargs):
        """
        Display raster data without projection
        figname: the file name to export map
        figsize: the size of map
        dpi: The resolution in dots per inch
        vmin and vmax define the data range that the colormap covers
        figname=None, figsize=None, dpi=300, vmin=None, vmax=None, 
                cax=True, dem_array=None, relocate=False, scale_ratio=1
        """
        fig, ax = gs.mapshow(raster_obj=self, **kwargs)
        return fig, ax
    
    def rankshow(self, **kwargs):
        """ Display water depth map in a range defined by (d_min, d_max)
        """
        fig, ax = gs.rankshow(self, **kwargs)
        return fig, ax
    
    def hillshade(self, **kwargs):
        """ Draw a hillshade map
        """
        fig, ax = gs.hillshade(self, **kwargs)
        return fig, ax

    def vectorshow(self, obj_y, **kwargs):
        """
        plot velocity map of U and V, whose values stored in two raster
        objects seperately
        """
        fig, ax = gs.vectorshow(self, obj_y, **kwargs)
        return fig, ax
#%%=========================== private functions ==============================
    def __osgeo2raster(self, obj_ds):
        """
        convert an osgeo dataset to a raster object
        """
        array = obj_ds.ReadAsArray()
        geo_trans_v = obj_ds.GetGeoTransform()
        projection = obj_ds.GetProjection()
        left = geo_trans_v[0]
        top = geo_trans_v[3]
        cellsize = geo_trans_v[1]
        nrows = obj_ds.RasterYSize
        ncols = obj_ds.RasterXSize
        xllcorner = left
        yllcorner = top - cellsize*nrows
        NODATA_value = obj_ds.GetRasterBand(1).GetNoDataValue()
        if NODATA_value is None:
            NODATA_value = -9999
        header = {'ncols':ncols, 'nrows':nrows,
                  'xllcorner':xllcorner, 'yllcorner':yllcorner,                  
                  'cellsize':cellsize, 'NODATA_value':NODATA_value}
        obj_new = Raster(array=array, header=header, projection=projection)
        return obj_new

    def __set_wkt_projection(self, epsg_code):
        """
        get coordinate reference system (crs) as Well Known Text (WKT) 
            from https://epsg.io
        epsg_code: the epsg code of a crs, e.g. BNG:27700, WGS84:4326
        return wkt text
        """
        import requests
        # access projection information
        wkt = requests.get('https://epsg.io/{0}.prettywkt/'.format(epsg_code))
        # remove spaces between charachters
        remove_spaces = wkt.text.replace(" ", "")
        # place all the text on one line
        output = remove_spaces.replace("\n", "")
        self.projection = output
        return output
    
#%%
def merge(obj_origin, obj_target, resample_method='bilinear'):
    """Merge the obj_origin to obj_target
    assign grid values in the origin Raster to the cooresponding grid cells in
    the target object. If cellsize are not equal, the origin Raster will be
    firstly resampled to the target object.
    obj_origin, obj_target: Raster objects
    """
    if obj_origin.header['cellsize'] != obj_target.header['cellsize']:
        obj_origin = obj_origin.resample(obj_target.header['cellsize'], 
                                   method=resample_method)
#    else:
#        obj_origin = self
    grid_x, grid_y = obj_origin.GetXYcoordinate()
    rows, cols = sp.map2sub(grid_x, grid_y, obj_target.header)
    ind_r = np.logical_and(rows >= 0, rows <= obj_target.header['nrows']-1)
    ind_c = np.logical_and(cols >= 0, cols <= obj_target.header['ncols']-1)
    ind = np.logical_and(ind_r, ind_c)
    ind = np.logical_and(ind, ~np.isnan(obj_origin.array))
    obj_output = copy.deepcopy(obj_target)
    obj_output.array[rows[ind], cols[ind]] = obj_origin.array[ind]
    return obj_output

def main():
    """Main function
    """
    print('Class to deal with raster data')

if __name__=='__main__':
    main()
    
    

