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

def fit_nl(x,a,b): 
     return b*(x**a)


def scaling(pop,param,plot=False):
    [a,err_a,b,err_b,r2,r,s] = [None]*7
    xy = np.array([[a,d] for a,d in zip(pop,param)])    
    if len(xy) == 0:
        return [a,err_a,b,err_b,r2,r,s]
    xy = xy[~np.isnan(xy).any(axis=1)]
    if len(xy) == 0:
        return [a,err_a,b,err_b,r2,r,s]
    xy = xy[np.logical_and(np.logical_and(xy[:,0] > 0, xy[:,1] > 0),np.sum(xy,axis=1)<np.inf)]
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

def scaling_onefixed(pop,param,plot=False):
    [a,err_a,b,err_b,r2,r,s] = [None]*7
    xy = np.array([[a,d] for a,d in zip(pop,param)])    
    if len(xy) == 0:
        return [a,err_a,b,err_b,r2,r,s]
    xy = xy[~np.isnan(xy).any(axis=1)]
    if len(xy) == 0:
        return [a,err_a,b,err_b,r2,r,s]
    xy = xy[np.logical_and(np.logical_and(xy[:,0] > 0, xy[:,1] > 0),np.sum(xy,axis=1)<np.inf)]
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

def scaling_nl(pop,param,plot=False):
    [a,err_a,b,err_b,r2,r,s] = [None]*7
    xy = np.array([[a,d] for a,d in zip(pop,param)])    
    if len(xy) == 0:
        return [a,err_a,b,err_b,r2,r,s]
    xy = xy[~np.isnan(xy).any(axis=1)]
    if len(xy) == 0:
        return [a,err_a,b,err_b,r2,r,s]
    xy = xy[np.logical_and(np.logical_and(xy[:,0] > 0, xy[:,1] > 0),np.sum(xy,axis=1)<np.inf)]
    if len(xy) <= 2:
        return [a,err_a,b,err_b,r2,r,s]
    popt,pcov=curve_fit(fit,np.log(xy[:,0]),np.log(xy[:,1]),p0=(1.0,0.0))
    a=popt[0] 
    b=popt[1]
    print([a,np.exp(b)])
    del popt, pcov
    xy = xy.astype(np.float128)    
    popt,pcov=curve_fit(fit_nl,xy[:,0],xy[:,1],p0=(a,np.exp(b))) 
    a=popt[0]; 
    err_a=np.sqrt(pcov[0,0])
    b=popt[1]
    err_b=np.sqrt(pcov[1,1])
    if np.isnan(err_a): err_a = 999
    if np.isnan(err_b): err_b = 999
    r2 = r2_score(xy[:,1],(xy[:,0]**a)*b)
    r,p = pearsonr((xy[:,0]**a)*b,xy[:,1])
    s,p = spearmanr((xy[:,0]**a)*b,xy[:,1])
    out = [a,err_a,b,err_b,r2,r,s]
    print(out)
    return out 


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
                        frac_houses = df_msa['bupl'].values/df_msa['all_bupl'].values
                        
                        adj_pop = pop*frac_houses
                    elif pop_calc == 'bufasum':
                        frac_bufasum = df_msa['bufasum'].values/df_msa['all_bufasum'].values
                        adj_pop = pop*frac_bufasum
                    elif pop_calc == 'bui':
                        frac_bui = df_msa['bui'].values/df_msa['all_bui'].values
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
                    #print([years,df_msa[param].values])
                    yearparams = df_msa[['year',param]].values#np.array([[y,p] for y,p in zip(years,df_msa[param].values)])
                    yearparams = yearparams[np.abs(yearparams[:,1])<np.inf]
                    years = yearparams[:,0]
                    params = yearparams[:,1]
                    a,err_a,r2 = trends(years,params)
                    trend_stats[param+'_a'].append(a)
                    trend_stats[param+'_err_a'].append(err_a)
                    trend_stats[param+'_r2'].append(r2)
            pd.DataFrame(data=trend_stats).to_csv(out_file,index=False)


def cross_analysis(df,patch,cumul,longitudinal,geo_completeness,temp_completeness,combined,pop_calc='buildings',non_linear=False):
    file_name = 'CrossSectionalScaling_TempComplete='+str(temp_completeness)+'_GeoComplete='+str(geo_completeness)+'_non_linear='+str(non_linear)+'_pop='+pop_calc+'.csv'
    if combined:
        file_name = 'CrossSectionalScaling_TempComplete='+str(temp_completeness)+'_GeoComplete='+str(geo_completeness)+'_pop_calc='+pop_calc+'_non_linear='+str(non_linear)+'_combined_patches.csv'

    if not longitudinal:
        df_complete = df.loc[(df['GeoComplete']>geo_completeness)&(df['TempComplete']>temp_completeness),]
        if cumul:
            years = list(range(1900,2020,10))+[2015]
            # patch vs city statistics
            out_file = 'CityStats/'+file_name
            if patch:
                out_file = 'CityPatchStats/'+file_name
            print(out_file)
            # population scaling laws:
            cross_scaling_stats = {'year':[]}
            for param in ['area','bui','dist','bufasum']:
                cross_scaling_stats[param+'_a']=[]
                cross_scaling_stats[param+'_err_a']=[]
                cross_scaling_stats[param+'_err_b']=[]
                cross_scaling_stats[param+'_b']=[]
                cross_scaling_stats[param+'_r2']=[]
                #cross_scaling_stats[param+'_r']=[]
                #cross_scaling_stats[param+'_s']=[]

            for year in years:
                # we only take data before 2015
                df_year = df_complete.loc[df_complete['year']==year,]
                # population is adjusted by fraction of houses in the patches
                if not combined:
                    pop = df_year['pop'].values
                    if pop_calc=='buildings':
                        frac_houses = df_year['patch_bupl'].values/df_year['all_bupl'].values
                        adj_pop = pop*frac_houses
                    elif pop_calc == 'bufasum':
                        frac_bufasum = df_year['bufasum'].values/df_year['all_bufasum'].values
                        adj_pop = pop*frac_bufasum
                    elif pop_calc == 'bui':
                        frac_bui = df_year['BUIs'].values/df_year['all_bui'].values
                        adj_pop = pop*frac_bui
                else:
                    if pop_calc=='buildings':
                        adj_pop = df_year['adj_pop_bupl'].values
                    elif pop_calc == 'bufasum':
                        adj_pop = df_year['adj_pop_bufa'].values
                    elif pop_calc == 'bui':
                        adj_pop = df_year['adj_pop_bui'].values

                # parameters of interest
                area = df_year['area'].values
                bui = df_year['BUIs'].values
                dist = df_year['distance'].values
                # calculate longitudinal scaling
                footprint = df_year['bufasum'].values

                if non_linear:
                    a_area,err_a_area,b_area,err_b_area,r2_area,r,s = scaling_nl(adj_pop,area)
                    a_bui,err_a_bui,b_bui,err_b_bui,r2_bui,r,s = scaling_nl(adj_pop,bui)
                    a_dist,err_a_dist,b_dist,err_b_dist,r2_dist,r,s = scaling_nl(adj_pop,dist)
                    a_fp,err_a_fp,b_fp,err_b_fp,r2_fp,r,s = scaling_nl(adj_pop,footprint)
                else:
                    print(scaling(adj_pop,area))
                    print(scaling(adj_pop,bui))
                    print(scaling(adj_pop,dist))
                    print(scaling(adj_pop,footprint))
                    a_area,err_a_area,b_area,err_b_area,r2_area,r,s = scaling(adj_pop,area)
                    a_bui,err_a_bui,b_bui,err_b_bui,r2_bui,r,s = scaling(adj_pop,bui)
                    a_dist,err_a_dist,b_dist,err_b_dist,r2_dist,r,s = scaling(adj_pop,dist)
                    a_fp,err_a_fp,b_fp,err_b_fp,r2_fp,r,s = scaling(adj_pop,footprint)


                cross_scaling_stats['year'].append(year)

                cross_scaling_stats['area_a'].append(a_area)
                cross_scaling_stats['area_err_a'].append(err_a_area)
                cross_scaling_stats['area_b'].append(b_area)
                cross_scaling_stats['area_err_b'].append(err_b_area)
                cross_scaling_stats['area_r2'].append(r2_area)
                #cross_scaling_stats['area_r'].append(r_area)
                #cross_scaling_stats['area_s'].append(s_area)

                cross_scaling_stats['bui_a'].append(a_bui)
                cross_scaling_stats['bui_err_a'].append(err_a_bui)
                cross_scaling_stats['bui_b'].append(b_bui)
                cross_scaling_stats['bui_err_b'].append(err_b_bui)
                cross_scaling_stats['bui_r2'].append(r2_bui)
                #cross_scaling_stats['bui_r'].append(r_bui)
                #cross_scaling_stats['bui_s'].append(s_bui)

                cross_scaling_stats['bufasum_a'].append(a_fp)
                cross_scaling_stats['bufasum_err_a'].append(err_a_fp)
                cross_scaling_stats['bufasum_b'].append(b_fp)
                cross_scaling_stats['bufasum_err_b'].append(err_b_fp)
                cross_scaling_stats['bufasum_r2'].append(r2_fp)
                #cross_scaling_stats['bufasum_r'].append(r_fp)
                #cross_scaling_stats['bufasum_s'].append(s_fp)

                cross_scaling_stats['dist_a'].append(a_dist)
                cross_scaling_stats['dist_err_a'].append(err_a_dist)
                cross_scaling_stats['dist_b'].append(b_dist)
                cross_scaling_stats['dist_err_b'].append(err_b_dist)
                cross_scaling_stats['dist_r2'].append(r2_dist)
                #cross_scaling_stats['dist_r'].append(r_dist)
                #cross_scaling_stats['dist_s'].append(s_dist)

            pd.DataFrame(data=cross_scaling_stats).to_csv(out_file,index=False)


geo_completeness,temp_completeness = [40,60]
cumul=True

for pop in ['buildings','bufasum','bui']:
    for patch in [True]:#,False]:
        print('patch: '+str(patch))
        if patch:
            directory = 'CityPatchStats/'
        else:
            directory = 'CityStats/'
        #for cumul in [True]:#[True,False]:
        for non_linear in [True]:#[False,True]:
            for combined in [True,False]:
                print('cumul: '+str(cumul))
                for longitudinal in [False]:#,True]:
                    stat_file = directory+'AllNetworkStats_cumulative='+str(cumul)+'.csv'
                    if combined and not longitudinal:
                        stat_file = directory+'AllNetworkStats_cumulative='+str(cumul)+'_allpatches.csv'
                    
                    df = pd.read_csv(stat_file)
                    #if combined and not longitudinal:
                    #    df['patch_bupl'] = df['bupl']
                    #    patch_stats = [ast.literal_eval(",".join(str(row['combined_patches']).replace('...','').replace('\n','').replace('\t','').replace('\r','').split()).replace('][','],[')) for ii,row in df.iterrows()]
                    #    msaids = [list(np.unique([pp[0]]) for p in patch_stats]
                    if not combined:
                        years = list(range(1900,2020,10))
                        all_list = []
                        for year in years:
                            all_dat = pd.read_csv('msa2val_'+str(year)+'.csv')
                            all_dat['year'] = [year]*len(all_dat)
                            all_list.append(all_dat)
                        all_bufa_bui = pd.concat(all_list).rename(columns={"bui": "all_bui", "bufa": "all_bufasum"})
                        bufa_data = pd.read_csv('MSA_BUFA_STATS.csv')
                        bufa_data['msaid'] = bufa_data['GEOID']
                        bufa_data['year'] = bufa_data['YEAR']
                        bufa_data['all_bufasum'] = bufa_data['BUFA_MSA_SUM']
                        all_bufa_bui = all_bufa_bui.drop(columns=['all_bufasum'])
                        all_bufa_bui = pd.merge(all_bufa_bui, bufa_data,on=['year','msaid'])
                        df = pd.merge(df,all_bufa_bui,on=['year','msaid'],how='outer')

                    if patch and not (combined and not longitudinal):
                        footprint_area_file = '../road_netw_multitemp_msa_level/MSBF_comparison/MSA_PATCH_BUFA_STATS.csv'
                        footprint_area_file_2015 = '../road_netw_multitemp_msa_level/MSBF_comparison/patches_2015_w_msbf_area_sum.csv'
                        footprint_area = pd.read_csv(footprint_area_file)
                        footprint_area['msaid'] = [val.split('_')[0] for val in footprint_area['uid'].values]
                        footprint_area['id'] = [val.split('_')[-1] for val in footprint_area['uid'].values]
                        footprint_area['bufasum'] = footprint_area['BUFA_PATCH_SUM']
                        for col in ['msaid','id','bufasum','YEAR']:
                            footprint_area[col] = footprint_area[col].astype(np.int64)
                        df = pd.merge(left=df,right=footprint_area,left_on=['msaid','id','year'],right_on=['msaid','id','YEAR'],how='inner')
                    elif not (combined and not longitudinal):
                        footprint_area_file = '../road_netw_multitemp_msa_level/MSBF_comparison/msbf_bufa_mutemp_cbsasum_all_lu.csv'
                        footprint_area = pd.read_csv(footprint_area_file)
                        df = pd.merge(left=df,right=footprint_area,left_on=['msaid','year'],right_on=['cbsa','year'],how='inner')
                    if non_linear and combined:
                        df = df.loc[df['adj_pop_bupl']<4*10**7,]
                    if not combined and longitudinal:
                        longitudinal_analysis(df,patch,cumul,longitudinal,geo_completeness,temp_completeness, combined,pop)
                    
                    cross_analysis(df,patch,cumul,longitudinal,geo_completeness,temp_completeness, combined,pop,non_linear)
