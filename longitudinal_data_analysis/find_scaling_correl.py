import pandas as  pd
import pickle as pk
import numpy as np
import geopandas as gpd
from geopy.distance import geodesic

cumul = True
print('CUMULATIVE = ',cumul)
#BUFA_cbsa = pd.read_csv('MSA_CBSA_BUFA_STATS.csv')
footprint_area_file = 'msbf_bufa_mutemp_cbsasum_all_lu.csv'
footprint_area = pd.read_csv(footprint_area_file)

df=pd.read_csv('CityStats/AllNetworkStats_cumulative='+str(cumul)+'.csv')#True.csv')
df = pd.merge(left=df,right=footprint_area,left_on=['msaid','year'],right_on=['cbsa','year'],how='left')

df_trend = pd.read_csv('CityStats/LongitudinalScaling.csv')
df['SA Type'] = df['MSA'].replace(0, 'uSA').replace(1,'MSA')
df['Region'] = df['region'].replace(1, 'NE').replace(2,'Midwest').replace(3,'South').replace(4,'West')
msa2pop = {msa:pop for msa,pop in df.loc[df['year']==2015,['msaid','pop']].values}
years = np.unique(df['year'].values)
for col in ['GeoComplete','TempComplete','SA Type','Region']:
    df_trend[col] = [df.loc[df['msaid']==row['msaid'],col].iloc[0] for n,row in df_trend.iterrows()]


for temp_complete,geo_coverage in [[0,0]]:#[[0,0],[80,80]]:
    df2=df.loc[(df['TempComplete']>temp_complete)&(df['GeoComplete']>geo_coverage),]
    df_trend_temp=df_trend.loc[(df_trend['TempComplete']>temp_complete)&(df_trend['GeoComplete']>geo_coverage),]
    df_trend_temp['pop'] = [msa2pop[msa] for msa in df_trend_temp['msaid'].values]
    
    # add msaid to name conversion
    PopPerCounty = pd.read_csv('state_pops/all_county_census_MSA_full.csv')
    select_cities = PopPerCounty['CBSA Code'].drop_duplicates().values
    CBSA_Names=PopPerCounty[['CBSA Code','CBSA Title']].dropna().drop_duplicates()
    df2['msa_name'] = [CBSA_Names.loc[CBSA_Names['CBSA Code']==city_id,'CBSA Title'].values[0] for city_id in df2['msaid'].values]
    
    ### specify relevant columns for t-SNE transform:
    #'area','dist','BUIs','BUFA','adj_pop'
    df2['adj_pop'] = df2['pop'].values*df2['patch_bupl'].values/df2['all_bupl'].values
    df_trend_temp['adj_pop'] = [df2.loc[df2['msaid']==row['msaid'],col].iloc[0] for n,row in df_trend_temp.iterrows()]
    relcols=['bui','bufasum','area','dist','adj_pop']

    #df2['SA Type'] = df2['MSA'].replace(0, 'uSA').replace(1,'MSA')
    #df2['Region'] = df2['region'].replace(1, 'NE').replace(2,'Midwest').replace(3,'South').replace(4,'West')
    #df2['num_nodes_per_area']=df2['num_nodes'].values/df2['area'].values*10**6
    #df2['num_edges_per_area']=df2['num_edges'].values/df2['area'].values*10**6
    #df2['distance_per_area']=df2['distance'].values/df2['area'].values*10**6
    #df2=df2[['msaid','msa_name','year','MSA','SA Type','Region']+relcols]
    ##### 
    # find patterns over time
    for msa in [-1]:#,0,1]:
        #df2 = df.loc[df['year']==years[yy],]
        df3 = df_trend_temp.copy(deep=True)
        if msa >= 0:
            df3 = df3.loc[(df3['MSA']==msa),]
        coord_file = 'msa_centroids_MSA='+str(msa)+'.pkl'
        msa_coordinates=pk.load(open(coord_file,'rb'))
        msaid_df3 = df3.groupby('msaid')
        for col in relcols:#[relcols[14]]:
            correl_file = col+'_spatial_correl_msa='+str(msa)+'_temp='+str(temp_complete)+'_geo='+str(geo_coverage)+'_scaling.csv'
            data_correl={'distance':[],'origin':[],'destination':[]}
            # find central coordinates for each msa
            val_coordinates=[]

            for n,line in df3.iterrows():
                msaid=line['msaid']
                if col != 'adj_pop':
                    diff = msaid_df3.get_group(msaid)[col+'_a'].values
                else:
                    diff = msaid_df3.get_group(msaid)[col].values
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
                        #shortest_len = min([len(val),len(val2)])
                        # assumptions: 
                        # - all data is complete after min_year
                        # - all np.nan is for the first or last few years (there is no gap/missing decade)
                        valid_val2 = val2[0]
                        #valid_val2 = valid_val2[~np.isnan(val2[-shortest_len:]) & ~np.isnan(val[-shortest_len:])]
                        valid_val = val[0]#[-shortest_len:]
                        #valid_val = valid_val[~np.isnan(val2[-shortest_len:]) & ~np.isnan(val[-shortest_len:])]
                        distance = geodesic(tuple(coord),tuple(coord2)).km
                        data_correl['distance'].append(distance)
                        data_correl['origin'].append(valid_val)
                        data_correl['destination'].append(valid_val2)
            data_correl = pd.DataFrame(data_correl)
            data_correl.to_csv(correl_file,index=False)
