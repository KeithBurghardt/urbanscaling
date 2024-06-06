import geopandas as gpd
import pickle as pk
import pandas as pd

df=pd.read_csv('CityStats/AllNetworkStats_cumulative=False.csv')
df['SA Type'] = df['MSA'].replace(0, 'uSA').replace(1,'MSA')


data = gpd.read_file('../RoadNetworks/outdata_spatial_merged/patches_allmsa_merged.shp')

wgs84_proj4string = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'  ###definition of lat lon coordinate system (wgs84)
data = data.to_crs(wgs84_proj4string)
msaid_data = data.groupby('msaid')
latest_year = 2015
data = data.loc[data['year']==latest_year,]
print(data.columns)
msaids = data['msaid'].drop_duplicates().values
for msa in [0,1,-1]:
    if msa>=0:
        msaids = df.loc[df['SA Type']==['uSA','MSA'][msa],'msaid'].drop_duplicates().values
    merged_shp = {msaid:msaid_data.get_group(msaid).unary_union.centroid for msaid in msaids}
    pk.dump(merged_shp,'msa_centroids_MSA='+str(msa)+'.pkl')











