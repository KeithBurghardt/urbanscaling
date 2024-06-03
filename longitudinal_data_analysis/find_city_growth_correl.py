import pandas as  pd
import pickle as pk
import numpy as np
import geopandas as gpd
from geopy.distance import geodesic

cumul = False
print('CUMULATIVE = ',cumul)
#BUFA_cbsa = pd.read_csv('MSA_CBSA_BUFA_STATS.csv')
footprint_area_file = 'msbf_bufa_mutemp_cbsasum_all_lu.csv'
footprint_area = pd.read_csv(footprint_area_file)

df=pd.read_csv('CityStats/AllNetworkStats_cumulative='+str(cumul)+'.csv')#True.csv')
df = pd.merge(left=df,right=footprint_area,left_on=['msaid','year'],right_on=['cbsa','year'],how='left')
temp_complete = 60
geo_coverage = 40
for temp_complete,geo_coverage in [[80,80],[0,0],[40,60]]:
    df=df.loc[(df['TempComplete']>temp_complete)&(df['GeoComplete']>geo_coverage),]
    # add msaid to name conversion
    PopPerCounty = pd.read_csv('state_pops/all_county_census_MSA_full.csv')
    select_cities = PopPerCounty.loc[PopPerCounty['POPESTIMATE2015']>3*10**5,'CBSA Code'].drop_duplicates().values
    CBSA_Names=PopPerCounty[['CBSA Code','CBSA Title']].dropna().drop_duplicates()
    df['msa_name'] = [CBSA_Names.loc[CBSA_Names['CBSA Code']==city_id,'CBSA Title'].values[0] for city_id in df['msaid'].values]

    ### specify relevant columns for t-SNE transform:
    #'area','distance','BUIs','BUFA','adj_pop'
    df['adj_pop'] = df['pop'].values*df['patch_bupl'].values/df['all_bupl'].values

    relcols=['area','distance','BUIs','bufasum','adj_pop','num_nodes_per_area','num_edges_per_area','distance_per_area','k_mean','k1','k4plus','entropy','mean_local_gridness','mean_local_gridness_max']
    years= list(range(1900,2020,10))+[2015]
    for feature in []:#'bufasum','area','distance','BUIs']:
        diff_feature = []
        for f,cbsa,y in df[[feature,'msaid','year']].values:
            if y == 1900:
                diff_feature.append(np.nan)
                continue
            year_pos = years.index(y)
            prev_year_f = df.loc[(df['msaid']==cbsa) & (df['year']==years[year_pos-1]),feature]
            if len(prev_year_f) == 0:
                diff_feature.append(np.nan)
                continue
            diff = f - prev_year_f.values[0]
            diff_feature.append(diff)
        df['diff_'+feature] = diff_feature
    relcols+=['pop']

    df['SA Type'] = df['MSA'].replace(0, 'uSA').replace(1,'MSA')
    df['Region'] = df['region'].replace(1, 'NE').replace(2,'Midwest').replace(3,'South').replace(4,'West')
    df['num_nodes_per_area']=df['num_nodes'].values/df['area'].values*10**6
    df['num_edges_per_area']=df['num_edges'].values/df['area'].values*10**6
    df['distance_per_area']=df['distance'].values/df['area'].values*10**6
    df=df[['msaid','msa_name','year','MSA','SA Type','Region']+relcols]
    ##### 
    # find patterns over time
    for msa in [-1]:#,0,1]:
        #df2 = df.loc[df['year']==years[yy],]
        df3 = df.copy(deep=True)#loc[df['year']==years[yy],]
        if msa >= 0:
            df3 = df3.loc[(df3['MSA']==msa),]
        coord_file = 'msa_centroids_MSA='+str(msa)+'.pkl'
        msa_coordinates=pk.load(open(coord_file,'rb'))
        msaid_df3 = df3.groupby('msaid')
        for col in relcols:#[relcols[14]]:
            correl_file = col+'_spatial_correl_msa='+str(msa)+'_temp='+str(temp_complete)+'_geo='+str(geo_coverage)+'_allyears_diff.csv'
            data_correl={'distance':[],'origin':[],'destination':[]}
            # find central coordinates for each msa
            val_coordinates=[]

            for n,line in df3.iterrows():
                msaid=line['msaid']
                diff = msaid_df3.get_group(msaid)[[col,'year']].sort_values(by='year')
                diff = diff[col].values
                if msaid in msa_coordinates.keys():
                    val_coordinates.append([diff,msa_coordinates[msaid]])
            dist_vals=[]
            for ii,vc in enumerate(val_coordinates):
                if ii % 100 == 0:
                    print(round(ii/len(val_coordinates)*100,2))
                for jj,vc2 in enumerate(val_coordinates):
                    if ii>jj:
                        val2,coord2 = vc2
                        coord2 = list(coord2.values[0].coords)[0][::-1]
                        val,coord = vc
                        coord = list(coord.values[0].coords)[0][::-1]
                        shortest_len = min([len(val),len(val2)])
                        # assumptions: 
                        # - all data is complete after min_year
                        # - all np.nan is for the first or last few years (there is no gap/missing decade)
                        valid_val2 = val2[-shortest_len:]
                        valid_val2 = valid_val2[~np.isnan(val2[-shortest_len:]) & ~np.isnan(val[-shortest_len:])]
                        valid_val = val[-shortest_len:]
                        valid_val = valid_val[~np.isnan(val2[-shortest_len:]) & ~np.isnan(val[-shortest_len:])]
                        distance = geodesic(tuple(coord),tuple(coord2)).km
                        data_correl['distance'].append(distance)
                        data_correl['origin'].append(valid_val)
                        data_correl['destination'].append(valid_val2)
            data_correl = pd.DataFrame(data_correl)
            data_correl.to_csv(correl_file,index=False)
