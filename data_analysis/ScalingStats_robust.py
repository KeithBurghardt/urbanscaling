import os,time,ast
import networkx as nx
import pandas as pd
import geopandas as gpd
import numpy as np
import geopandas as gp
from shapely.geometry import Point
import fiona
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from osgeo import gdal
import inspect
#import osmnx
from scipy.special import entr
from scipy.optimize import curve_fit
from scipy import stats
from scipy.stats import norm
from scipy.stats import spearmanr
from sklearn.neighbors import KernelDensity
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from sklearn.metrics import r2_score

def fit(x,a,b): 
     return a*x+b
def scaling(pop,param,plot=False):
    [a,err_a,b,err_b,r2,r,s] = [None]*7
    xy = np.array([[a,d] for a,d in zip(pop,param)])    
    if len(xy) == 0:
        return [a,err_a,b,err_b,r2,r,s]
    xy = xy[~np.isnan(xy).any(axis=1)]
    if len(xy) == 0:
        return [a,err_a,b,err_b,r2,r,s]
    xy = xy[np.logical_and(xy[:,0] > 0, xy[:,1] > 0)]
    if len(xy) <= 2:
        return [a,err_a,b,err_b,r2,r,s]

    popt,pcov=curve_fit(fit,np.log(xy[:,0]),np.log(xy[:,1]),p0=(1.0,0.0)) 
    a=popt[0]; 
    err_a=np.sqrt(pcov[0,0])
    b=popt[1]
    err_b=np.sqrt(pcov[1,1])
    r2 = r2_score(np.log(xy[:,1]),a*np.log(xy[:,0])+b)
    r,p = pearsonr(a*np.log(xy[:,0])+b,np.log(xy[:,1]))
    s,p = spearmanr(a*np.log(xy[:,0])+b,np.log(xy[:,1]))
    return [a,err_a,b,err_b,r2,r,s] 

def trends(years,param):
    [a,err_a,r2] = [None]*3
    xy = np.array([[a,d] for a,d in zip(years,param)])
    if len(xy) == 0:
        return [a,err_a,r2]
    xy = xy[~np.isnan(xy).any(axis=1)]
    if len(xy) <= 2:
        return [a,err_a,r2]
    popt,pcov=curve_fit(fit,xy[:,0],xy[:,1],p0=(1.0,0.0)) 
    a=popt[0]; 
    b=popt[1]
    err_a=np.sqrt(pcov[0,0])
    r2 = r2_score(a*xy[:,0]+b,xy[:,1])
    return [a,err_a,r2] 



def longitudinal_analysis(df,patch,cumul,longitudinal,geo_completeness,temp_completeness, combined,pop_calc='buildings'):
    out_file = 'CityStats/LongitudinalScaling_TempComplete='+str(temp_completeness)+'_GeoComplete='+str(geo_completeness)+'_pop_calc='+pop_calc+'.csv'
                    
    df_complete = df.loc[(df['GeoComplete']>geo_completeness)&(df['TempComplete']>temp_completeness),]
    if longitudinal and not patch:
        if cumul:
            for start_year,end_year in [[1900,2015],[1900,1960],[1960,2015]]:
                out_file = 'CityStats/LongitudinalScaling_TempComplete='+str(temp_completeness)+'_GeoComplete='+str(geo_completeness)+'_pop_calc='+pop_calc+'_start='+str(start_year)+'_end='+str(end_year)+'.csv'
                # population scaling laws:
                long_scaling_stats = {'msaid':[]}
                long_scaling_stats['GeoComplete'] = []
                long_scaling_stats['TempComplete'] = []
                for param in ['area','bui','dist','bufasum']:
                    long_scaling_stats[param+'_a']=[]
                    long_scaling_stats[param+'_err_a']=[]
                    long_scaling_stats[param+'_b']=[]
                    long_scaling_stats[param+'_err_b']=[]
                    long_scaling_stats[param+'_r2']=[]
                for msaid in df_complete['msaid'].drop_duplicates().values:
                    # we only take data before 2015
                    df_msa = df_complete.loc[(df_complete['msaid']==msaid)&(df_complete['year']<end_year) & (df_complete['year']>= start_year),]
                    if len(df_msa) == 0: continue
                    long_scaling_stats['GeoComplete'].append(df_msa['GeoComplete'].values[0])
                    long_scaling_stats['TempComplete'].append(df_msa['TempComplete'].values[0])
                    # population is adjusted by fraction of houses in the patches
                    pop = df_msa['pop'].values
                    if pop_calc=='buildings':
                        frac_houses = df_msa['patch_bupl'].values/df_msa['all_bupl'].values
                        adj_pop = pop*frac_houses
                    elif pop_calc == 'bufasum':
                        frac_bufasum = df_msa['patch_bufasum'].values/df_msa['all_bufasum'].values
                        adj_pop = pop*frac_bufasum
                    elif pop_calc == 'bui':
                        frac_bui = df_msa['patch_bui'].values/df_msa['all_bui'].values
                        adj_pop = pop*frac_bui

                    # parameters of interest
                    area = df_msa['area'].values
                    bui = df_msa['BUIs'].values
                    dist = df_msa['distance'].values
                    footprint = df_msa['bufasum'].values
                    # calculate longitudinal scaling
                    a_area,err_a_area,b_area,err_b_area,r2_area,r,s = scaling(adj_pop,area)
                    a_bui,err_a_bui,b_bui,err_b_bui,r2_bui,r,s = scaling(adj_pop,bui)
                    a_dist,err_a_dist,b_dist,err_b_dist,r2_dist,r,s = scaling(adj_pop,dist)
                    a_fp,err_a_fp,b_fp,err_b_fp,r2_fp,r,s = scaling(adj_pop,footprint)

                    long_scaling_stats['msaid'].append(msaid)
                    long_scaling_stats['area_a'].append(a_area)
                    long_scaling_stats['area_err_a'].append(err_a_area)
                    long_scaling_stats['area_b'].append(b_area)
                    long_scaling_stats['area_err_b'].append(err_b_area)
                    long_scaling_stats['area_r2'].append(r2_area)
                    long_scaling_stats['bui_a'].append(a_bui)
                    long_scaling_stats['bui_err_a'].append(err_a_bui)
                    long_scaling_stats['bui_b'].append(b_bui)
                    long_scaling_stats['bui_err_b'].append(err_b_bui)
                    long_scaling_stats['bui_r2'].append(r2_bui)
                    long_scaling_stats['dist_a'].append(a_dist)
                    long_scaling_stats['dist_err_a'].append(err_a_dist)
                    long_scaling_stats['dist_b'].append(b_dist)
                    long_scaling_stats['dist_err_b'].append(err_b_dist)
                    long_scaling_stats['dist_r2'].append(r2_dist)
                    long_scaling_stats['bufasum_a'].append(a_fp)
                    long_scaling_stats['bufasum_err_a'].append(err_a_fp)
                    long_scaling_stats['bufasum_b'].append(b_fp)
                    long_scaling_stats['bufasum_err_b'].append(err_b_fp)
                    long_scaling_stats['bufasum_r2'].append(r2_fp)
                pd.DataFrame(data=long_scaling_stats).to_csv(out_file,index=False)


        else:
            # calculate whether yearly network statistics are increasing or decreasing
            out_file = 'CityStats/CityStatTrends.csv'
            trend_stats = {'msaid':[]}
            for param in ['num_nodes_per_area','num_edges_per_area','distance_per_area','k_mean','k1','k4plus','entropy','mean_local_gridness','mean_local_gridness_max']:
                trend_stats[param+'_a']=[]
                trend_stats[param+'_err_a']=[]
                trend_stats[param+'_r2']=[]
            for msaid in df_complete['msaid'].drop_duplicates().values:
                # we only take data before 2015 
                df_msa = df_complete.loc[(df_complete['msaid']==msaid)&(df_complete['year']<2015),]
                # assume areas are collected sequentially
                areas = df_msa['area'].values
                d_areas = areas[1:]-areas[:-1]
                # remove 1900 data
                df_msa = df_msa.iloc[1:] #df.loc[(df['msaid']==msaid)&(1900<df['year'])&(df['year']<2015),]
                years = df_msa['year'].values
                num_nodes_per_area = df_msa['num_nodes'].values/d_areas*10**6
                num_edges_per_area = df_msa['num_edges'].values/d_areas*10**6
                distance_per_area = df_msa['distance'].values/d_areas*10**6
                df_msa.loc[:,'num_nodes_per_area']=num_nodes_per_area[:]
                df_msa.loc[:,'num_edges_per_area']=num_edges_per_area[:]
                df_msa.loc[:,'distance_per_area']=distance_per_area[:]
                #df_msa = df_msa.dropna()
                trend_stats['msaid'].append(msaid)
                for param in ['num_nodes_per_area','num_edges_per_area','distance_per_area','k_mean','k1','k4plus','entropy','mean_local_gridness','mean_local_gridness_max']:
                    yearparams = df_msa[['year',param]].values#np.array([[y,p] for y,p in zip(years,df_msa[param].values)])
                    yearparams = yearparams[np.abs(yearparams[:,1])<np.inf]
                    years = yearparams[:,0]
                    params = yearparams[:,1]
                    a,err_a,r2 = trends(years,params)
                    trend_stats[param+'_a'].append(a)
                    trend_stats[param+'_err_a'].append(err_a)
                    trend_stats[param+'_r2'].append(r2)
            pd.DataFrame(data=trend_stats).to_csv(out_file,index=False)


def cross_analysis(df,patch,cumul,longitudinal,geo_completeness,temp_completeness,combined,file_name,pop_calc='buildings'):

    if not longitudinal:
        df_complete = df.loc[(df['GeoComplete']>geo_completeness)&(df['TempComplete']>temp_completeness),]
        if cumul:
            years = list(range(1900,2020,10))+[2015]
            # patch vs city statistics
            out_file = 'CityStats/'+file_name
            if patch:
                out_file = 'CityPatchStats/'+file_name
            # population scaling laws:
            cross_scaling_stats = {'year':[]}
            for param in ['area','bui','dist','bufasum']:
                cross_scaling_stats[param+'_a']=[]
                cross_scaling_stats[param+'_err_a']=[]
                cross_scaling_stats[param+'_b']=[]
                cross_scaling_stats[param+'_err_b']=[]
                cross_scaling_stats[param+'_r2']=[]
                cross_scaling_stats[param+'_r']=[]
                cross_scaling_stats[param+'_s']=[]

            for year in years:
                # we only take data before 2015
                df_year = df_complete.loc[df_complete['year']==year,]
                pop = df_year['pop'].values
                if pop_calc=='buildings':
                    frac_houses = df_year['patch_bupl'].values/df_year['all_bupl'].values
                    adj_pop = pop*frac_houses
                elif pop_calc == 'bufasum':
                    frac_bufasum = df_year['patch_bufasum'].values/df_year['all_bufasum'].values  
                    adj_pop = pop*frac_bufasum
                elif pop_calc == 'bui':
                    frac_bui = df_year['patch_bui'].values/df_year['all_bui'].values  
                    adj_pop = pop*frac_bui

                # parameters of interest
                area = df_year['area'].values
                bui = df_year['BUIs'].values
                dist = df_year['distance'].values
                footprint = df_year['bufasum'].values
                #plt.plot(adj_pop,footprint,'k.')
                #plt.yscale('log')
                #plt.xscale('log')
                #plt.show()

                # calculate longitudinal scaling
                a_area,err_a_area,b_area,err_b_area,r2_area,r_area,s_area = scaling(adj_pop,area)
                a_fp,err_a_fp,b_fp,err_b_fp,r2_fp,r_fp,s_fp = scaling(adj_pop,footprint)
                #a_bui,err_a_bui,b_bui,err_b_bui,r2_bui,r_bui,s_bui = [None]*5
                a_bui,err_a_bui,b_bui,err_b_bui,r2_bui,r_bui,s_bui = scaling(adj_pop,bui)
                a_dist,err_a_dist,b_dist,err_b_dist,r2_dist,r_dist,s_dist = scaling(adj_pop,dist)
                cross_scaling_stats['year'].append(year)
                cross_scaling_stats['area_a'].append(a_area)
                cross_scaling_stats['area_err_a'].append(err_a_area)
                cross_scaling_stats['area_b'].append(b_area)
                cross_scaling_stats['area_err_b'].append(err_b_area)
                cross_scaling_stats['area_r2'].append(r2_area)
                cross_scaling_stats['area_r'].append(r_area)
                cross_scaling_stats['area_s'].append(s_area)
                cross_scaling_stats['bui_a'].append(a_bui)
                cross_scaling_stats['bui_err_a'].append(err_a_bui)
                cross_scaling_stats['bui_b'].append(b_bui)
                cross_scaling_stats['bui_err_b'].append(err_b_bui)
                cross_scaling_stats['bui_r2'].append(r2_bui)
                cross_scaling_stats['bui_r'].append(r_bui)
                cross_scaling_stats['bui_s'].append(s_bui)                
                cross_scaling_stats['bufasum_a'].append(a_fp)
                cross_scaling_stats['bufasum_err_a'].append(err_a_fp)
                cross_scaling_stats['bufasum_b'].append(b_fp)
                cross_scaling_stats['bufasum_err_b'].append(err_b_fp)
                cross_scaling_stats['bufasum_r2'].append(r2_fp)
                cross_scaling_stats['bufasum_r'].append(r_fp)
                cross_scaling_stats['bufasum_s'].append(s_fp)
                cross_scaling_stats['dist_a'].append(a_dist)
                cross_scaling_stats['dist_err_a'].append(err_a_dist)
                cross_scaling_stats['dist_b'].append(b_dist)
                cross_scaling_stats['dist_err_b'].append(err_b_dist)
                cross_scaling_stats['dist_r2'].append(r2_dist)
                cross_scaling_stats['dist_r'].append(r_dist)
                cross_scaling_stats['dist_s'].append(s_dist)
            pd.DataFrame(data=cross_scaling_stats).to_csv(out_file,index=False)

pop = 'buildings'#
#pop = 'bufasum'
#pop = 'bui'

radii=[500,1000,2000]  
thresholds=[0.05,0.10,0.20]
patch_percentiles=[0.8,0.9,0.95] 
directory='robustness_checks/'

cumul=True
longitudinal = False
combined=False

geo_completeness,temp_completeness=[40,60]
#  minimum distance between clusters - I used 1 cellsize (250m) and 2 cellsizes (500m)
for min_dist in [1,2]:
    df_list = []
    for year in [1900,1960,2010]:
        df_temp = pd.read_csv(directory+'cca_sum_df'+str(min_dist)+'_'+str(year)+'.csv')
        df_temp['year'] = [year]*len(df_temp)
        df_list.append(df_temp)
    
    df = pd.concat(df_list)
    print(df.columns)
    # remove cells < 1km
    df = df.loc[df['clustersize_cells']>=16,]
    #cca_cluster_id	msaid	bupr	bupl	bua	bui	bufa	kmroad
    df['id'] = df['cca_cluster_id']
    df['patch_bupr'] = df['bupr']
    df['patch_bupl'] = df['bupl']        
    df['bufasum'] = df['bufa']
    df['BUIs'] = df['bui']
    df['area'] = df['bua']
    df['distance'] = df['kmroad']
    df['msaid'] = [str(msa).zfill(5) for msa in df['msaid'].astype(str).values]
    print(len(df))
    print(len(df.dropna()))

    all_bupl = pd.read_csv('CityPatchStats/AllNetworkStats_cumulative=True.csv')[['msaid','year','all_bupl','all_bupr','TempComplete','GeoComplete','pop']]
    all_bupl['msaid'] = [str(msa).zfill(5) for msa in all_bupl['msaid'].astype(str).values]
    cca_df= pd.merge(df,all_bupl,on = ['msaid','year'])
    patch=True
    print('patch: '+str(patch))
    file_name = 'CrossSectionalScaling_TempComplete='+str(temp_completeness)+'_GeoComplete='+str(geo_completeness)+'_min_dist='+str(min_dist)+'_cca.csv'
    cross_analysis(cca_df,patch,cumul,longitudinal,geo_completeness,temp_completeness, combined,file_name,pop)        


for geo_completeness,temp_completeness in [[0,0],[40,60],[0,0],[80,80]]:
    #if False:
    print([geo_completeness,temp_completeness])
    for r in radii:
        print(r)
        for t in thresholds:
            print(t)
            file = directory+'segment_stats_sensitivity_'+str(r)+'_'+str(t)+'.csv'
            df = pd.read_csv(file)
            df['id'] = df['Id']
            df['patch_bupr'] = df['bupr_sum']
            df['patch_bupl'] = df['bupl_sum']
            
            df['bufasum'] = df['bufa_sum']
            df['BUIs'] = df['bui_sum']
            df['area'] = df['bua_sum']
            df['distance'] = df['roadkm_sum']
            footprint_area_file = '../road_netw_multitemp_msa_level/MSBF_comparison/MSA_PATCH_BUFA_STATS.csv'
            footprint_area = pd.read_csv(footprint_area_file)
            footprint_area['msaid'] = [val.split('_')[0] for val in footprint_area['uid'].values]
            footprint_area['id'] = [val.split('_')[-1] for val in footprint_area['uid'].values]
            #footprint_area['bufasum'] = footprint_area['BUFA_PATCH_SUM']
            for col in ['msaid','id','YEAR']:
                footprint_area[col] = footprint_area[col].astype(np.int64)
            print(df.columns)
            print(footprint_area.columns)
            df = pd.merge(left=df,right=footprint_area,left_on=['msaid','id','year'],right_on=['msaid','id','YEAR'],how='inner')
            
            for pp in patch_percentiles:
                print(pp)
                pp_data = df.loc[df['percentile_selection']==pp,]
                all_bupl = pd.read_csv('CityPatchStats/AllNetworkStats_cumulative=True.csv')[['id','msaid','year','all_bupl','all_bupr','TempComplete','GeoComplete','pop']]
                pp_data= pd.merge(pp_data,all_bupl,left_on = ['id','msaid','year'],right_on=['id','msaid','year'])
                patch=True
                print('patch: '+str(patch)) 
                file_name = 'CrossSectionalScaling_TempComplete='+str(temp_completeness)+'_GeoComplete='+str(geo_completeness)+'_rad='+str(r)+'_thresh='+str(t)+'_patch_percent='+str(pp)+'_robust.csv'

                cross_analysis(pp_data,patch,cumul,longitudinal,geo_completeness,temp_completeness, combined,file_name,pop)
