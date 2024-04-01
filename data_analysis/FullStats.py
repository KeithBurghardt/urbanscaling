import os,time,ast
import networkx as nx
import pandas as pd
import numpy as np
import geopandas as gp
from shapely.geometry import Point
import fiona
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import gdal
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


directory='/Volumes/Keith Network Hard Drive/RoadNetworks/outdata_roads/'
# sub-regions to add later
# what region is city?
city_regions = pd.read_csv('network_stats/all_county_census_MSA.csv')[['CBSA Code','REGION']].drop_duplicates().values
city_regions={city_id:region for city_id,region in city_regions}

# is city MSA or uSA
city_MSA = pd.read_csv('network_stats/all_county_census_MSA.csv')[['CBSA Code','Metropolitan/Micropolitan Statistical Area']].drop_duplicates().values
city_MSA ={city_id:int(MSA=='Metropolitan Statistical Area') for city_id,MSA in city_MSA}


# completeness
# cols: 'GEOID', temp_completeness,geographic_coverage	
completeness=pd.read_csv('network_stats/MSA_UNC_STATS.csv')

#MSAID ->  BUI, AREA 
# GEOID -> BUI_MSA_SUM
df_areaBUI = pd.read_csv('network_stats/segment_stats_ALL_incl_bui.csv')
#msaid,Id,year -> area, bui_sum



# combine all city stats into one file
bin_width = 2*np.pi/60;

## 'patch_bupr':[],'patch_bupl':[],'patch_bua':[]: find patch_* for all IDs that become new ID, find patch_* of new ID, then subtract
## 'area': 
## FullStats.py: # if not cumul: find IDs of previous years that become current ID, find previous area, find current, subtract
## nx.read_file
## bigger goal: find what IDs become what?? use networkx nodes to find IDs with (complete) overlap. Area from BUI, etc.
def find_past_patches(directory,years,city_id):
    nodes_per_year = {'year':[],'id':[],'nodes':[]}
    for y_index,year in enumerate(years):
        year_col = 'CENSUS'+str(year)+'POP'
        if year == 2015:
            year_col = 'POPESTIMATE2015'
        city_file = directory+'roads_'+str(city_id)+'_1000_005_'+str(year)+'.shp'
        old_year = year -10
        if year == 2015:
            old_year = year-5
        old_city_file = directory+'roads_'+str(city_id)+'_1000_005_'+str(old_year)+'.shp'
        if os.path.exists(city_file):# and os.path.exists(old_city_file):
            data=gp.read_file(city_file)
            for ID in data['Id'].drop_duplicates().values.flatten():
                in_file = directory+'roads_'+str(city_id)+'_1000_005_'+str(year)+'_id='+str(ID)+'.shp'
                if not os.path.exists(in_file):
                    data.loc[data['Id']==ID,].to_file(filename= in_file)
                G = nx.read_shp(in_file).to_undirected()
                nodes=set([n for n in G.nodes])
                nodes_per_year['year'].append(year)
                nodes_per_year['id'].append(ID)
                nodes_per_year['nodes'].append(set(nodes))
    nodes_per_year = pd.DataFrame(nodes_per_year)
    past_decade_IDs = []
    for n,line in nodes_per_year.iterrows():
        ID = line['id']
        year = line['year']
        past_year = year-10
        # no past IDs
        #if past_year < 1900:
        #    past_decade_IDs.append([])
        #    continue
        # if we are at 2015
        if past_year == 2005:
            past_year = 2010
        current_nodes = line['nodes']
        parent_ids=[]
        past_df = nodes_per_year.loc[nodes_per_year['year']==past_year,]
        if len(past_df) > 0:
            previous_ids = {line['id']:line['nodes'] for n,line in past_df.iterrows()}
            for key in previous_ids.keys():
                intersect = previous_ids[key].intersection(current_nodes)
                if len(intersect) > 1: # more than 1 of nodes are the same
                    #if len(intersect) < len(previous_ids[key]):
                    parent_ids.append(key)
        past_decade_IDs.append(parent_ids)
    nodes_per_year['prev_ids']=past_decade_IDs
    nodes_per_year = nodes_per_year.drop(columns=['nodes'])
    return nodes_per_year

def create_yearly_stats_from_cumul(df,file,directory):
    new_file = file.replace('stats_new=True','stats_new=False')#.replace('_no2degree=True','')
     # use an older file if the newer one doesnt exist
    if not os.path.exists(new_file):
        new_file = file.replace('stats_new=True','stats_new=False').replace('_no2degree=True','')
    cumul_df = pd.read_csv(new_file)
    # if this is true, then the non-cumulative data is bad
    #print(df.columns)
    years = np.unique(df['year'].sort_values().values)
    # make sure statistics are cumulative: look at cumulative df
    #print(df['msaid'].values[0])        
    if patch:
        past_ids = find_past_patches(directory,years,city_id)
        #print(past_ids['prev_ids'].values[0])

    for col in ['patch_bupr','patch_bupl','patch_bua','area','BUIs','distance','all_bupr','all_bupl','all_bua']:
        yearly_stats = []
        for n,row in df.iterrows():
            cid = row['msaid']
            year = row['year']
            past_year = year - 10
            # correct for 2015 year
            if past_year == 2005:
                past_year = 2010
            past_stat = 0
            if patch:
                ID = row['id'] 
                prev_ids = past_ids.loc[(past_ids['id']==ID)&(past_ids['year']==year),'prev_ids'].values[0]
                if col in cumul_df.columns:    
                    current_stat = cumul_df.loc[(cumul_df['msaid']==cid)&(cumul_df['id']==ID)&(cumul_df['year']==year),col].values[0]
                    if len(prev_ids) > 0:
                        past_stat = np.sum([cumul_df.loc[(cumul_df['msaid']==cid)&(cumul_df['id']==pid)&(cumul_df['year']==past_year),col].values for pid in prev_ids])
                else:#elif col not in ['all_bupr','all_bupl','all_bua']:
                    current_stat = df.loc[(df['msaid']==cid)&(df['id']==ID)&(df['year']==year),col].values[0]
                    if len(prev_ids) > 0:
                        past_stat = np.sum([df.loc[(df['msaid']==cid)&(df['id']==pid)&(df['year']==past_year),col].values for pid in prev_ids])
                #else:
                #    current_stat = cumul_df.loc[(cumul_df['msaid']==cid)&(df['id']==ID)&(cumul_df['year']==year),col].values[0]
                #    past_candidate_stat = cumul_df.loc[(cumul_df['msaid']==cid)&(df['id']==ID)&(cumul_df['year']==past_year),col].values[0]
                if col in ['all_bupr','all_bupl','all_bua']:
                    current_stat = cumul_df.loc[(cumul_df['msaid']==cid)&(cumul_df['id']==ID)&(cumul_df['year']==year),col].values[0]
                    past_stat = 0

            else:
                if col in cumul_df.columns:    
                    current_stat = cumul_df.loc[(cumul_df['msaid']==cid)&(cumul_df['year']==year),col].values[0]
                    past_candidate_stat = cumul_df.loc[(cumul_df['msaid']==cid)&(cumul_df['year']==past_year),col].values
                else:
                    current_stat = df.loc[(df['msaid']==cid)&(df['year']==year),col].values[0]
                    past_candidate_stat = df.loc[(df['msaid']==cid)&(df['year']==past_year),col].values  
                if len(past_candidate_stat) > 0:
                    past_stat = past_candidate_stat[0]

                if col in ['all_bupr','all_bupl','all_bua']:
                    current_stat = cumul_df.loc[(cumul_df['msaid']==cid)&(cumul_df['year']==year),col].values[0]
                    past_stat = 0

            yearly_stats.append(current_stat-past_stat)
        df[col] = yearly_stats
    return df

second_century=True
for patch in [False]:#[False,True]:
    print('patch: '+str(patch))
    if patch:
        cumul_files=sorted(list(glob('CityPatchStats/NetworkStatistics_*stats_new=False*_no2degree*')))
        yearly_files=sorted(list(glob('CityPatchStats/NetworkStatistics_*stats_new=True*_no2degree*')))
    else:
        if second_century:
            cumul_files= list(glob('CityStats/NetworkStatistics_city-id=*_stats_new=True_no2degree=True_second_century=True.csv'))
        else:
            # first add BUI, Area
            cumul_files = list(glob('CityStats/NetworkStatistics_*stats_new=False*_no2degree*'))
            cumul_files = [f for f in cumul_files if 'second_century=True' not in f]
            yearly_files = list(glob('CityStats/NetworkStatistics_*stats_new=True*_no2degree*'))
            yearly_files = [f for f in yearly_files if 'second_century=True' not in f]

    for cumul in [True]:#,True]:
        print('cumul: '+str(cumul))
        if cumul:
            files = cumul_files[:]
        else:
            files = yearly_files[:]
        
        out_file = 'CityStats/AllNetworkStats_cumulative='+str(cumul)+'_second-century='+str(second_century)+'.csv'
        if patch:
            out_file = 'CityPatchStats/AllNetworkStats_cumulative='+str(cumul)+'.csv'
        all_df_list = []
        files = sorted(files)
        for nn,file in enumerate(files):
            if nn % 5 == 0:
                print(round(nn/len(files)*100,2))
            df = pd.read_csv(file)
            if len(df) == 0:
                continue
            city_id = df['msaid'].values[0]
            if len(df) == 0: continue
            if patch:
                areas = [df_areaBUI.loc[(df_areaBUI['msaid']==row['msaid'])&(df_areaBUI['Id']==row['id'])&(df_areaBUI['year']==row['year']),'area'].values[0] for n,row in df.iterrows()]
                # if not cumul: find IDs of previous years that become current ID, find previous area, find current, subtract
                BUIs = [df_areaBUI.loc[(df_areaBUI['msaid']==row['msaid'])&(df_areaBUI['Id']==row['id'])&(df_areaBUI['year']==row['year']),'bui_sum'].values[0] for n,row in df.iterrows()]
                push = 4
            else:
                areas = [np.sum(df_areaBUI.loc[(df_areaBUI['msaid']==row['msaid'])&(df_areaBUI['year']==row['year']),'area'].values) for n,row in df.iterrows()]
                BUIs = [np.sum(df_areaBUI.loc[(df_areaBUI['msaid']==row['msaid'])&(df_areaBUI['year']==row['year']),'bui_sum'].values) for n,row in df.iterrows()]
                push = 3
            df['area'] = areas
            df['BUIs'] = BUIs

            # replace 'patch_bupr','patch_bupl','patch_bua','all_bupr','all_bupl','all_bua','distance'
            if not cumul:
                # create cumulative statistics for 'patch_bupr','patch_bupl','patch_bua','all_bupr','all_bupl','all_bua','distance'
                # ignore the decade-by-decade statistics previously recorded
                df = create_yearly_stats_from_cumul(df,file,directory)
            MSAs = [city_MSA[city_id]]*len(df)
            df['MSA'] = MSAs

            Regions = [city_regions[city_id]]*len(df)
            TempComplete = [completeness.loc[completeness['GEOID']==city_id,'temp_completeness'].values[0]]*len(df)
            GeoComplete = [completeness.loc[completeness['GEOID']==city_id,'geographic_coverage'].values[0]]*len(df)
            df['region'] = Regions
            df['TempComplete']=TempComplete
            df['GeoComplete']=GeoComplete
            for remove_col in ['Unnamed: 0','index']:
                if remove_col in df.columns:
                    df = df.drop(remove_col,axis='columns')
            cols = df.columns.tolist()
            cols = cols[:push]+cols[-6:]+cols[push:-6]
            df = df[cols]
            #print(df)
            all_df_list.append(df)
        all_df = pd.concat(all_df_list)
        all_df.to_csv(out_file,index=False)
