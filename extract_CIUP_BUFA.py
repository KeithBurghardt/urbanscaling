import pandas as pd
import numpy as np

combined_patches = pd.read_csv('../road_netw_multitemp_msa_level/allpatches_merged_cbsa.csv')
print(len(combined_patches.loc[combined_patches['year']==2015,]))
def fit(x,a,b):
     return a*x+b

outfile = 'CityPatchStats/AllNetworkStats_cumulative=True_allpatches.csv'

if not os.path.exists(outfile):
    df = pd.read_csv('CityPatchStats/AllNetworkStats_cumulative=True.csv')
    footprint_area_file = '../road_netw_multitemp_msa_level/MSBF_comparison/MSA_PATCH_BUFA_STATS.csv'
    footprint_area_file_2015 = '../road_netw_multitemp_msa_level/MSBF_comparison/patches_2015_w_msbf_area_sum.csv'
    footprint_area = pd.read_csv(footprint_area_file)
    footprint_area['msaid'] = [val.split('_')[0] for val in footprint_area['uid'].values]
    footprint_area['id'] = [val.split('_')[-1] for val in footprint_area['uid'].values]
    footprint_area['bufasum'] = footprint_area['BUFA_PATCH_SUM']
    footprint_area_2015 = pd.read_csv(footprint_area_file_2015)
    footprint_area_2015['YEAR'] = [2015]*len(footprint_area_2015)
    footprint_area_2015['msaid'] = [val.split('_')[0] for val in footprint_area_2015['uid2'].values]
    footprint_area_2015['id'] = [val.split('_')[-1] for val in footprint_area_2015['uid2'].values]
    footprint_area_2015['bufasum'] = footprint_area_2015['SUM']
    footprint_area = pd.concat([footprint_area,footprint_area_2015])

    for col in ['msaid','id','bufasum','YEAR']:
        footprint_area[col] = footprint_area[col].astype(np.int64)
    df = pd.merge(left=df,right=footprint_area,left_on=['msaid','id','year'],right_on=['msaid','id','YEAR'],how='inner')

    merged_patches = [[p.split('_') for p in ast.literal_eval(patches.replace('\n','').replace(' ',','))] for patches in combined_patches['uid2_of_dissolved_patches'].values]
    combined_stats = {'year':[],'adj_pop':[],'area':[],'BUIs':[],'bufasum':[],'distance':[],'GeoComplete':[],'TempComplete':[],'region':[],'combined_patches':[]}
    for ii,patches in enumerate(merged_patches):
        all_patch_stats = {'year':0,'adj_pop':0,'area':0,'BUIs':0,'bufasum':0,'distance':0,'GeoComplete':0,'TempComplete':0}
        patches = np.array(patches).astype(int)
        #statistics = df.loc[(df['year']==year),]
        #print(df.columns)
        #print(statistics.columns)
        patches = np.unique(patches, axis=0)
        valid_patches = 0
        region = ''
        for patch in patches:
            msaid,year,pid = patch
            statistics = df.loc[(df['year']==year),]
            statistics = statistics.loc[(statistics['msaid']==msaid) & (statistics['id'] == pid),['year','all_bupl','patch_bupl','pop','area','BUIs','bufasum','distance','GeoComplete','Temp$
            statistics = statistics.loc[(statistics['pop']>0)&(statistics['patch_bupl']>0),]

            if len(statistics) > 0:
                region = statistics['region'].values[0]
                valid_patches+=1
                adj_pop = statistics['patch_bupl'].values/statistics['all_bupl'].values*statistics['pop'].values[0]
                statistics['adj_pop'] = adj_pop
                for key in ['adj_pop','area','BUIs','bufasum','distance','GeoComplete','TempComplete']:
                    all_patch_stats[key] += statistics[key].values[0]
        #print(statistics.columns)
        if all_patch_stats['adj_pop'] > 0:
            #print(combined_stats.keys())
            all_patch_stats['year'] = year
            # take the average
            all_patch_stats['GeoComplete'] /= valid_patches
            all_patch_stats['TempComplete'] /= valid_patches
            combined_stats['combined_patches'].append(patches)
            combined_stats['region'].append(region)
            for key in combined_stats.keys():
                if key not in ['region','combined_patches']:
                    combined_stats[key].append(all_patch_stats[key])

    combined_stats = pd.DataFrame(combined_stats)
    combined_stats.to_csv(outfile)

