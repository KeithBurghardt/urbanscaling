# -*- coding: utf-8 -*-
"""
Created on Sat Dec 23 22:48:44 2023

@author: Johannes H. Uhl, University of Colorado Boulder, USA.
"""
## SENSITIVITY ANALYSIS OF GENERALIZED BUILT-UP AREAS
## creates generalized BU areas (GBUAs) for a range of parameters, as shapefiles
## creates a shapefile of patches with attributes per parameter combination
## params are: 
    ## thresholds
    ## patch_percentiles
    ## radii
## outputs: a csv file per scenario.

## call the script with arcpy (propy):
## "C:\Program Files\ArcGIS\Pro\bin\Python\Scripts\propy" <script name>
## note that an ArcGIS Pro license is required

import arcpy
import time
import os
import numpy as np
from arcpy import env
from arcpy.sa import *
from scipy.stats import rankdata
arcpy.CheckOutExtension("Spatial")
arcpy.env.compression='LZW'
arcpy.env.overwriteOutput=True
import pandas as pd

## gridded input data from HISDAC-US etc:
fbuy='./IN_DATA/00_FBUY_COMPOSITE_semidecadal_V2.tif' #HISDAC-US, https://doi.org/10.7910/DVN/PKJ90M, 
bupr_dummy='./IN_DATA/BUPR/BUPR_XXXX.tif' #HISDAC-US, https://doi.org/10.7910/DVN/YSWMDR, 
bupl_dummy='./IN_DATA/BUPL/BUPL_XXXX.tif' #HISDAC-US, https://doi.org/10.7910/DVN/SJ213V,
bua_dummy='./IN_DATA/BUA/BUA_XXXX.tif' #HISDAC-US, https://doi.org/10.7910/DVN/J6CYUJ
bui_dummy='./IN_DATA/BUI/BUI_XXXX.tif' #HISDAC-US, https://doi.org/10.7910/DVN/1WB9E4
bufa_dummy='./IN_DATA/hist_msbf_bufa_accum_XXXX.tif'  #HISDAC-US, https://doi.org/10.7910/DVN/HXQWNJ,
roadkm_dummy='./IN_DATA/gridcell_stats_kmroad_1km_all_cbsas.tif' #https://doi.org/10.6084/m9.figshare.19593496.v1
## CBSA geometry and population data:
msa_shp='./IN_DATA/tl_2019_us_cbsa.shp' #https://catalog.data.gov/dataset/tiger-line-shapefile-2019-nation-u-s-current-metropolitan-statistical-area-micropolitan-statist
county_msa_pop_csv='./IN_DATA/all_county_census_MSA_full.csv'

tempfolder = './TEMP/'
tempgdb= './TEMP/tempdata.gdb' # to be created manually, empty file geodatabase
outfolder = './OUT_DATA/'

years = [1900,1950,2010]
radii=[500,1000,2000]  #in m 
tresholds=[0.05,0.10,0.20] #,
patch_percentiles=[0.8,0.9,0.95]

do_sensitivity = True

####################################################################

arcpy.env.outputCoordinateSystem = arcpy.Describe(fbuy).spatialReference
arcpy.env.workspace=outfolder
arcpy.env.scratchWorkspace=tempfolder

county_msa_pop_df = pd.read_csv(county_msa_pop_csv)
county_msa_pop_df['STATE']=county_msa_pop_df.FIPS.map(str).str.zfill(5).str.slice(0,2)

msas_all = [[row[0],row[1]] for row in arcpy.da.SearchCursor(msa_shp, ['GEOID','NAME'])]    
msas_all = sorted(msas_all)

if do_sensitivity:
    
    for radius_m in radii:
        for treshold in tresholds:
            patch_stats_df=pd.DataFrame()
            msacount=0
            curr_treshold_val = treshold*(np.pi*radius_m*radius_m/(250.0*250.0))
            for msa in msas_all:                 
                msastart = time.time()
                msacount+=1
                msa_geoid = msa[0]
                msa_name = msa[1]
                ### get states that MSA covers.
                states_involved = np.unique(county_msa_pop_df[county_msa_pop_df['CBSA Code']==int(msa_geoid)]['STATE'].values)
                ### get counties that MSA covers.            
                counties_involved = np.unique(county_msa_pop_df[county_msa_pop_df['CBSA Code']==int(msa_geoid)]['FIPS'].values)
                counties_involved = [str(int(x)).zfill(5) for x in counties_involved]	
                clipshp = tempfolder+os.sep+'temp_msa.shp'
                arcpy.Select_analysis(msa_shp,clipshp,""" "GEOID" = '%s' """%msa_geoid)               
                refyear=years[0]  
                init_minsizes=[]
                for year in years:
                    try:
                        bupr_in = bupr_dummy.replace('XXXX',str(year))
                        bupr_clipped = ExtractByMask(bupr_in, clipshp)                       
                    except:
                        continue
                    bupl_in = bupl_dummy.replace('XXXX',str(year))    
                    bua_in = bua_dummy.replace('XXXX',str(year))
                    bui_in = bui_dummy.replace('XXXX',str(year))
                    bufa_in = bufa_dummy.replace('XXXX',str(year))
                    roadkm_in = roadkm_dummy.replace('XXXX',str(year))
                    
                    bua_rast = Con(bupr_clipped>0,1,0)  
                    arcpy.CopyRaster_management(bua_rast,outfolder+os.sep+'bua_%s.tif' %year)  
                    
                    #### spatial generalization
                    neighborhood = NbrCircle(radius_m, "MAP")
                    outFocalStatistics = FocalStatistics(bua_rast, neighborhood, "SUM","")
                    
                    #### density threshold
                    outFocalStatistics_tresh = Con(outFocalStatistics>curr_treshold_val,1)
                    arcpy.CopyRaster_management(outFocalStatistics_tresh,outfolder+os.sep+'generalized_extent_%s_%s_%s.tif' %(int(100*treshold),radius_m,year))    
                    print(treshold,year)
                    outPolygons = outfolder+os.sep+'generalized_extent_poly_%s_%s_%s_%s.shp' %(int(100*treshold),radius_m,year,msa_geoid)
                    
                    #### segmentation of generalized surface
                    arcpy.RasterToPolygon_conversion(outFocalStatistics_tresh, outPolygons, "NO_SIMPLIFY", 'VALUE')        
                    arcpy.AddField_management(outPolygons,"area","Double")
                    expression1 = "{0}".format("!SHAPE.area@SQUAREMETERS!")         
                    arcpy.CalculateField_management(outPolygons, "area", expression1, "PYTHON", ) 

                    try:
                        ### get number of structures per segment
                        zonalstats_tbl_bupr = tempgdb+os.sep+'bupr_sums'
                        outZSaT = ZonalStatisticsAsTable(outPolygons, 'Id', bupr_in, zonalstats_tbl_bupr, "NODATA", "SUM")                
                    except:
                        ### in case of error, no settlements exist in year.
                        try:
                            refyear=years[years.index(year)+1]  ### this was missing 01-2021   
                        except:
                            continue
                        continue
                    try:
                        zonalstats_tbl_bupl = tempgdb+os.sep+'bupl_sums'
                        outZSaT = ZonalStatisticsAsTable(outPolygons, 'Id', bupl_in, zonalstats_tbl_bupl, "NODATA", "SUM")                
    
                        zonalstats_tbl_bui = tempgdb+os.sep+'bui_sums'
                        outZSaT = ZonalStatisticsAsTable(outPolygons, 'Id', bui_in, zonalstats_tbl_bui, "NODATA", "SUM")                
            
                        zonalstats_tbl_bua = tempgdb+os.sep+'bua_sums'
                        outZSaT = ZonalStatisticsAsTable(outPolygons, 'Id', bua_in, zonalstats_tbl_bua, "NODATA", "SUM")  
                        
                        zonalstats_tbl_bufa = tempgdb+os.sep+'bufa_sums'
                        outZSaT = ZonalStatisticsAsTable(outPolygons, 'Id', bufa_in, zonalstats_tbl_bufa, "NODATA", "SUM")                
                                    
                        zonalstats_tbl_roadkm = tempgdb+os.sep+'roadkm_sums'
                        outZSaT = ZonalStatisticsAsTable(outPolygons, 'Id', roadkm_in, zonalstats_tbl_roadkm, "NODATA", "SUM")  
                    except:
                        try:
                            refyear=years[years.index(year)+1]  ### this was missing 01-2021   
                        except:
                            continue
                        continue                        
                                    
                    #### read all into pandas and generate segment size statistics               
                    
                    poly_df = pd.DataFrame([[row[0],row[1]] for row in arcpy.da.SearchCursor(outPolygons, ['Id','area'])])
                    poly_df.columns=['Id','area']
                    
                    bupr_sums_df = pd.DataFrame([[row[0],row[1]] for row in arcpy.da.SearchCursor(zonalstats_tbl_bupr, ['Id','SUM'])])
                    bupr_sums_df.columns=['Id','bupr_sum']
                                    
                    bupl_sums_df = pd.DataFrame([[row[0],row[1]] for row in arcpy.da.SearchCursor(zonalstats_tbl_bupl, ['Id','SUM'])])
                    bupl_sums_df.columns=['Id','bupl_sum']

                    bui_sums_df = pd.DataFrame([[row[0],row[1]] for row in arcpy.da.SearchCursor(zonalstats_tbl_bui, ['Id','SUM'])])
                    bui_sums_df.columns=['Id','bui_sum']
                    
                    bua_sums_df = pd.DataFrame([[row[0],row[1]] for row in arcpy.da.SearchCursor(zonalstats_tbl_bua, ['Id','SUM'])])
                    bua_sums_df.columns=['Id','bua_sum']

                    bufa_sums_df = pd.DataFrame([[row[0],row[1]] for row in arcpy.da.SearchCursor(zonalstats_tbl_bufa, ['Id','SUM'])])
                    bufa_sums_df.columns=['Id','bufa_sum']
                    
                    roadkm_sums_df = pd.DataFrame([[row[0],row[1]] for row in arcpy.da.SearchCursor(zonalstats_tbl_roadkm, ['Id','SUM'])])
                    roadkm_sums_df.columns=['Id','roadkm_sum']
                                    
                    segment_df = poly_df.merge(bupr_sums_df,on='Id').merge(bupl_sums_df,on='Id').merge(bui_sums_df,on='Id').merge(bua_sums_df,on='Id').merge(bufa_sums_df,on='Id').merge(roadkm_sums_df,on='Id')
                    segment_df['msaid']=msa_geoid
                    segment_df['msaname']=msa_name
                    segment_df['year']=year
                    segment_df['radius_m']=radius_m
                    segment_df['treshold']=treshold

                    segment_df = segment_df[segment_df['bupr_sum']>0] ## remove uninhabited segments
                    
                    segment_df['bupl_rank']=rankdata(segment_df.bupl_sum.values,method='ordinal')
                    segment_df['bupl_pcntl']=segment_df.bupl_sum.rank(pct=True)

                    for patch_percentile in patch_percentiles:                                    
                        try:
                            if year==refyear: ## smallest segment in minyear is the min size in all subsequent years.
                                init_minsize = np.min(segment_df[segment_df.bupl_pcntl>patch_percentile].bupl_sum.values)
                                init_minsizes.append([patch_percentile,init_minsize])
                            print(init_minsize)
                        except:
                            try:
                                refyear=years[years.index(year)+1]  ### this was missing 01-2021   
                            except:
                                continue
                            continue                               
                        try:
                            ## minsize is never reached
                            init_minsize = init_minsizes[patch_percentiles.index(patch_percentile)][1]
                        except:
                            continue
                        ##### make a selection  which segments belong to GBUA. ###############
                        segment_df = segment_df[np.logical_or(segment_df.bupl_pcntl>patch_percentile,segment_df.bupl_sum > init_minsize)]                            
                        segment_df['percentile_selection']=patch_percentile
                        segment_df['init_minsize_bupl']=init_minsize
                        patch_stats_df = patch_stats_df.append(segment_df)               
                        print(msacount,'patch stats done',radius_m,treshold,patch_percentile,msa_geoid,msa_name,states_involved,year)
                        patch_stats_df.to_csv(tempfolder+os.sep+'segment_stats_sensitivity_temp_%s_%s.csv' %(radius_m,treshold),index=False)                    
                
            patch_stats_df.to_csv(tempfolder+os.sep+'segment_stats_sensitivity_%s_%s.csv' %(radius_m,treshold),index=False)                    
                                   
                                        
