import pandas as pd
import geopandas as gpd

# get all BUFA, BUI for MSA
msa_gdf = gpd.read_file('tl_2019_us_cbsa/tl_2019_us_cbsa.shp').to_crs(epsg=3857)

def rescale_vals(df):
    standard_area = 62500
    df['area'] = df['geometry'].area
    df['new_dn'] = df['DN'].values*df['area'].values/standard_area
    df['geometry'] =  df['geometry'].centroid    
    return df

years = list(range(1900,2020,10))+[2015]
for year in years:
    print(year)
    bufa = gpd.read_file('/Volumes/Keith_Burghardt_EHD/RoadNetworks/bufa/hisdac_us_bufa_'+str(year)+'.shp')
    print(len(bufa))
    print('bufa uploaded')
    bufa = rescale_vals(bufa).to_crs(epsg=3857)
    print(len(bufa))
    print('bufa rescaled')
    bui = gpd.read_file('/Volumes/Keith_Burghardt_EHD/RoadNetworks/bui/BUI/data/BUI_'+str(year)+'.shp')
    print('bui uploaded')
    bui = rescale_vals(bui).to_crs(epsg=3857)
    print('finished')
    bufa_within = gpd.sjoin(bufa, msa_gdf, how='inner', predicate='within')
    bui_within = gpd.sjoin(bui, msa_gdf, how='inner', predicate='within')
    bufa_msaid = set(bufa_within['GEOID'].drop_duplicates().values.tolist())    
    bui_msaid = set(bui_within['GEOID'].drop_duplicates().values.tolist())    
    bufa_within = bufa_within.groupby('GEOID')
    bui_within = bui_within.groupby('GEOID')
    print(bufa_msaid)
    print(bui_msaid)
    print(msa_gdf['GEOID'].values.tolist())
    msa2val={'msaid':[],'bui':[],'bufa':[]}    
    for ii,row in msa_gdf.iterrows():
        if ii % 10 == 0:
            print(round(ii/len(msa_gdf)*100,2))
        msa = row['geometry']
        msaid = row['GEOID']
        msa2val['msaid'].append(msaid)
        bufa_sum = 0
        if msaid in bufa_msaid:
            bufa_sum = bufa_within.get_group(msaid)['new_dn'].values.sum()
        bui_sum = 0
        if msaid in bui_msaid:        
            bui_sum = bui_within.get_group(msaid)['new_dn'].values.sum()
        msa2val['bufa'].append(bufa_sum)
        msa2val['bui'].append(bui_sum)
    msa2val = pd.DataFrame(msa2val)
    msa2val.to_csv('msa2val_'+str(year)+'.csv',index=False)
