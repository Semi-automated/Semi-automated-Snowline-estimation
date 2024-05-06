# -*- coding: utf-8 -*-
"""
@author: Parul
#Snowline estimation from Sentinel-2 and Landsat
# Before running the follwing code please replace 'input/output folder location' with your input/output folder locations, define total number of images and define buffer size according to the sensor (Sentinel-2/Landsat)
"""


import rasterio
import numpy as np
import geopandas as gd
import fiona
import glob
import rasterio.env
from rasterio.features import shapes
import matplotlib.pyplot as plt
import skimage.color
import skimage.filters
import cv2
import statistics as stats
import pandas as pd

#keep same extent for DEM, NIR band and SWIR band files
#reading DEM
DEM=rasterio.open('input DEM location with file name and extension (.tif)','r')
dem=DEM.read(1).astype('float64')

#reading all the images (NIR and SWIR bands)
NIRlist = glob.glob('input folder location'+ '/**/*NIR.tif', recursive=True)
SWIRlist = glob.glob('input folder location'+ '/**/*SWIR.tif', recursive=True)

thresholds=[]
SLA_values=[]
Max_SLA=[]
Min_SLA=[]
for i in range(number of images+1):
    b=NIRlist[i]
    B=SWIRlist[i]
    NIR1=rasterio.open(b)
    NIR2=NIR1.read(1).astype('float32')
    NIR2[(NIR2 == 65535) ] = np.nan
    figure_size=3
    NIR3=cv2.medianBlur(NIR2, figure_size)
    NIR=NIR3.astype('float64')
    SWIR1=rasterio.open(B)
    SWIR2=SWIR1.read(1).astype('float32')
    SWIR2[(SWIR2 == 65535) ] = np.nan
    SWIR3=cv2.medianBlur(SWIR2, figure_size)
    SWIR=SWIR3.astype('float64')
    Bandratio=NIR/SWIR
    print(Bandratio)
    Bandratio_image=rasterio.open('output folder/Bandratio_image'+str(i)+'.tif','w',driver='Gtiff',
                width=NIR1.width,height=NIR1.height,
                count=1,
                crs=NIR1.crs,
                transform=NIR1.transform,
                dtype='float64'
                )
    Bandratio_image.write(Bandratio,1)
    Bandratio_image.close()
    
    #scaling factor for both Landsat and Sentinel is 1000
    scaling_factor = 1000
    NIRB=NIR/scaling_factor
    
    # snow and ice
    snowice=np.where((Bandratio>1.5),1,0)
    snowice=snowice.astype('float64')#to convert from int to float
    
    # perform automatic thresholding
    th = skimage.filters.threshold_otsu(NIRB[np.isfinite(NIRB)])
    print("Found automatic threshold th = {}.".format(th))
    
    # create a histogram of the NIR band
    histogram, bin_edges = np.histogram(NIRB, bins=256, range=(0.0, 1.0))
    fig, ax = plt.subplots()
    plt.plot(bin_edges[0:-1], histogram)
    plt.title("Histogram")
    plt.xlabel("Reflectance")
    plt.ylabel("pixels")
    plt.xlim(0, 1.0)
    ax.axvline(x=th, color='r', linestyle='dashed')
    
    #snow on ice and bareland with others
    snow_on_ice=np.where((snowice==1) & (NIRB>th),
                     1,
                     0)
    snow_on_ice=snow_on_ice.astype('float64')#to convert from int to float
    
    #Exposed ice
    exposed_ice=snowice-snow_on_ice
    exposed_ice=exposed_ice.astype('float64')
   
    F_snow_on_iceImage=rasterio.open('output folder location/F_snow_on_iceImage'+str(i)+'.tif','w',driver='Gtiff',
                width=NIR1.width,height=NIR1.height,
                count=1,
                crs=NIR1.crs,
                transform=NIR1.transform,
                dtype='int32'
                )
    F_snow_on_iceImage.write(snow_on_ice,1)
    F_snow_on_iceImage.close()
   
    O_exposed_iceImage=rasterio.open('output folder location/O_exposed_iceImage'+str(i)+'.tif','w',driver='Gtiff',
                width=NIR1.width,height=NIR1.height,
                count=1,
                crs=NIR1.crs,
                transform=NIR1.transform,
                dtype='int32'
                )
    O_exposed_iceImage.write(exposed_ice,1)
    O_exposed_iceImage.close()
    
    #raster to polygon and buffer
    buffer_size=5 or 15 #for Sentinel-2 5m and for Landsat 15m
    with rasterio.open('output folder location/F_snow_on_iceImage'+str(i)+'.tif') as src1:
        image1 = src1.read()
        mask1 = image1 != 0
    results1 = ({'properties': {'raster_val': v}, 'geometry': s}
    for i, (s, v)in enumerate(shapes(image1,mask=mask1,transform=src1.transform)))
    with fiona.open('output folder location/snow'+str(i)+'.shp', 'w',driver='Shapefile',crs=src1.crs,schema={'properties': [('raster_val', 'int')],'geometry': 'Polygon'}) as dst1:
        dst1.writerecords(results1)
    Snow=gd.read_file('output folder location/snow'+str(i)+'.shp')
    Snow['geometry']=Snow['geometry'].buffer(buffer_size)
    Snow.to_file('output folder location/buffer_S'+str(i)+'.shp')
    
    # same for ice buffer
    with rasterio.open('output folder location/O_exposed_iceImage'+str(i)+'.tif') as src:
        image = src.read()
        mask = image != 0
    results = ({'properties': {'raster_val': v}, 'geometry': s}
    for i, (s, v) in enumerate(shapes(image,mask=mask,transform=src.transform)))
    with fiona.open('output folder location/ice'+str(i)+'.shp', 'w', driver='Shapefile',crs=src.crs,schema={'properties': [('raster_val', 'int')],'geometry': 'Polygon'}) as dst:
        dst.writerecords(results)
    Ice=gd.read_file('output folder location/ice'+str(i)+'.shp')
    Ice['geometry']=Ice['geometry'].buffer(buffer_size)
    Ice.to_file('output folder location/buffer_I'+str(i)+'.shp')
    
    #converting to geodataframe
    from geopandas import GeoDataFrame, overlay
    points = overlay(Snow, Ice, how='intersection', keep_geom_type=True)
    
    #intersection
    intersect = points.unary_union
    
    # extract DEM values
    from rasterio.mask import mask
    masked, mask_transform = mask(dataset=DEM, shapes=[intersect])
    masked=masked.astype('float64')
    masked[(masked == 65535) ] = np.nan 
    maskedImage=rasterio.open('output folder location/maskedImage'+str(i)+'.tif','w',driver='Gtiff',
                width=NIR1.width,height=NIR1.height,
                count=1,
                crs=NIR1.crs,
                transform=NIR1.transform,
                dtype='float64'
                )
    maskedImage.write(masked)
    maskedImage.close()
    print('initial shape: {}'.format(masked.shape))
    masked = masked.squeeze()
    print('final shape: {}'.format(masked.shape))
    
    #clipping to the glaciers boundaries -100m buffer
    masked2=rasterio.open('output folder location/maskedImage'+str(i)+'.tif')
    Glaciers=gd.read_file('input glacier boundaries shapefile location with extension (.shp)')
    Glacier_B=Glaciers.buffer(-100)
    Glaciers_MP = Glacier_B.unary_union
    masked_clip, mask_transform = mask(dataset=masked2, shapes=[Glaciers_MP])
    masked_clip=masked_clip.astype('float64')
    masked_clip[(masked_clip == 0) ] = np.nan 
    maskedImage_clip=rasterio.open('output folder location/maskedImage_clip'+str(i)+'.tif','w',driver='Gtiff',
                width=NIR1.width,height=NIR1.height,
                count=1,
                crs=NIR1.crs,
                transform=NIR1.transform,
                dtype='float64'
                )
    maskedImage_clip.write(masked_clip)
    maskedImage_clip.close()
    print('initial shape: {}'.format(masked_clip.shape))
    masked_clip = masked_clip.squeeze()
    print('final shape: {}'.format(masked_clip.shape))
    
    #for plotting box plot
    m=masked_clip.flatten()
    m = m[~np.isnan(m)]
    
    df=pd.DataFrame(m)
    figure=plt.boxplot(df)
    plt.ylabel('SLA Elevation (m a.s.l.)')
    from matplotlib.cbook import boxplot_stats
    stats=boxplot_stats(df.values)
    figure.keys()
    from operator import itemgetter
    q1 = list(map(itemgetter('q1'), stats))
    q3 = list(map(itemgetter('q3'), stats))
    whislo = list(map(itemgetter('whislo'), stats))
    whishi = list(map(itemgetter('whishi'), stats))
    masked_clip_th=np.where((masked_clip<q3[0])&(masked_clip>whislo[0]),masked_clip,np.nan)
    masked_clip_thImage=rasterio.open('output folder location/masked_clip_thImage'+str(i)+'.tif','w',driver='Gtiff',
                width=NIR1.width,height=NIR1.height,
                count=1,
                crs=NIR1.crs,
                transform=NIR1.transform,
                dtype='float64',
                nodata=65535
                )
    masked_clip_thImage.write(masked_clip_th,1)
    masked_clip_thImage.close()
    
    #removing isolated pixels with sieve filter, prepare binary image
    filtered=np.where((masked_clip_th>0)&(masked_clip_th != np.nan),1,0)
    filteredImage=rasterio.open('output folder location/filteredImage'+str(i)+'.tif','w',driver='Gtiff',
                width=NIR1.width,height=NIR1.height,
                count=1,
                crs=NIR1.crs,
                transform=NIR1.transform,
                dtype='float64',
                nodata=65535
                )
    filteredImage.write(filtered,1)
    filteredImage.close()
    from osgeo import gdal
    Image = gdal.Open('output folder location/filteredImage'+str(i)+'.tif', 1)  # open image in read-write mode
    Band = Image.GetRasterBand(1)
    gdal.SieveFilter(srcBand=Band, maskBand=None, dstBand=Band, threshold=50, connectedness=8, callback=gdal.TermProgress_nocb)
    del Image, Band  # close the datasets
    
    #again extract valid DEM values after sieve filter
    f=rasterio.open('output folder location/filteredImage'+str(i)+'.tif')
    f1=f.read(1)
    f1[(f1 == 0) ] = np.nan
    SLP_value=np.where((f1==1)&(dem != 65535),dem,np.nan)
    SLPImage=rasterio.open('output folder location/SLPImage'+str(i)+'.tif','w',driver='Gtiff',
                width=NIR1.width,height=NIR1.height,
                count=1,
                crs=NIR1.crs,
                transform=NIR1.transform,
                dtype='float64',
                nodata=65535
                )
    SLPImage.write(SLP_value,1)
    SLPImage.close()
    
    #SLA calculation
    SLA=np.nanmedian(SLP_value)
    max_SLA=np.nanmax(SLP_value)
    min_SLA=np.nanmin(SLP_value)
    print(SLA)
    print(max_SLA)
    print(min_SLA)
    SLA_values.append(SLA)
    Max_SLA.append(max_SLA)
    Min_SLA.append(min_SLA)
    thresholds.append(th)

dp = pd.DataFrame(list(zip(NIRlist, SLA_values, Max_SLA, Min_SLA, thresholds)),
               columns =['NIRlist', 'SLA_values', 'Max_SLA', 'Min_SLA', 'thresholds'])
with pd.ExcelWriter('output folder location/SLAs.xlsx') as writer:
	dp.to_excel(writer, sheet_name='SLAs')
