import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import spearmanr
import matplotlib.cm as cm
from scipy.optimize import curve_fit



plt.style.use('default')
def fit(x,a,b):
     return a*x+b
scaling_params = ['resident_pop','area','built_up_area','length_total']
sp=cm.jet
labelfonts = {'fontname':'Arial','fontsize':15}
markers=['^','o','s','D','*','v','<','>','1','2','3','4']
markers=['^','o','s','D','*','v','<','>']#,'1','2','3','4']
params = {"ytick.color" : "k",
          "xtick.color" : "k",
          "axes.labelcolor" : "k",
          "axes.edgecolor" : "k"}
plt.rcParams.update(params)
scaling_data={}
fig, axes = plt.subplots(1, 3,figsize=(17,5))
coeffs_subregion={'feature':[],'region':[],'a':[],'err_a':[],'b':[],'err_b':[],'n':[]}
agg_coeffs_subregion={'feature':[],'a':[],'err_a':[],'b':[],'err_b':[]}
for ff,feature in enumerate(scaling_params[1:]):
    
    plt_xy=boeing_data_ww.copy()
    popt,pcov=curve_fit(fit,np.log(plt_xy['resident_pop'].values),np.log(plt_xy[feature].values),p0=(0.0,0.0)) 
    a=popt[0]; 
    err_a=np.sqrt(pcov[0,0])
    #yr_stats=[year,a,err_a]
    b=popt[1]; 
    err_b=np.sqrt(pcov[1,1])
    agg_coeffs_subregion['a'].append(a)
    agg_coeffs_subregion['err_a'].append(err_a)
    agg_coeffs_subregion['b'].append(b)
    agg_coeffs_subregion['err_b'].append(err_b)
    agg_coeffs_subregion['feature'].append(feature)
    subregions = boeing_data_ww['world_subregion'].drop_duplicates().values
    regions = boeing_data_ww['world_region'].drop_duplicates().values
    list_regions = []
    for r in regions:
        list_regions+= list(boeing_data_ww.loc[boeing_data_ww['world_region']==r,'world_subregion'].drop_duplicates().sort_values().values)

    for ii,region in enumerate(reversed(list_regions)):
        marker=markers[ii % len(markers)]
        plt_xy=boeing_data_ww.loc[boeing_data_ww['world_subregion']==region,]
        popt,pcov=curve_fit(fit,np.log(plt_xy['resident_pop'].values),np.log(plt_xy[feature].values),p0=(0.0,0.0)) 
        a=popt[0]; 
        err_a=np.sqrt(pcov[0,0])
        #yr_stats=[year,a,err_a]
        b=popt[1]; 
        err_b=np.sqrt(pcov[1,1])
        coeffs_subregion['region'].append(region)
        coeffs_subregion['a'].append(a)
        coeffs_subregion['err_a'].append(err_a)
        coeffs_subregion['b'].append(b)
        coeffs_subregion['err_b'].append(err_b)
        coeffs_subregion['feature'].append(feature)
        coeffs_subregion['n'].append(len(plt_xy))
        perr = np.sqrt(np.diag(pcov))
        pop = np.array([1,10**9])
        axes[ff].plot(pop,np.exp(b)*pop**a,color=sp((ii-1)/18))
    for ii,region in enumerate(reversed(list_regions)):
        marker=markers[ii % len(markers)]
        plt_xy=boeing_data_ww.loc[boeing_data_ww['world_subregion']==region,]
        # plt_xy['resident_pop'],plt_xy[feature]
        axes[ff].scatter(plt_xy['resident_pop'],plt_xy[feature],color=sp((ii-1)/18),marker=marker,alpha=0.3,label=region)# vary color by world region/country
    axes[ff].set_yscale('log')
    axes[ff].set_xscale('log')
    #plt.legend(ncol=2, bbox_to_anchor=(1.0, 0.8), loc='upper left')
    axes[ff].set_xlim([10**4.5,10**8])
    axes[ff].set_xlabel('Population',**labelfonts)
    if feature == 'length_total':
        axes[ff].set_ylabel('Road Length ($km$)',**labelfonts)
        axes[ff].set_ylim([10**3,10**8.5])
        axes[ff].set_aspect(0.6)
    else:
        if feature == 'area':
            axes[ff].set_ylabel('Urban Center Area ($km^2$)',**labelfonts)
        if feature == 'built_up_area':
            axes[ff].set_ylabel('Builtup Area ($km^2$)',**labelfonts)

        axes[ff].set_ylim([10**-0.1,10**4.1])
        axes[ff].set_aspect(0.8)
    
    axes[ff].tick_params(axis='both', which='major', labelsize=14)
axes[0].legend(loc='upper right', bbox_to_anchor=(0.5, -0.5))
plt.tight_layout()
plt.savefig('WW_Scaling.pdf',transparent=True)
plt.show()
plt.close()




scaling_data={}
fig, axes = plt.subplots(1, 3,figsize=(17,5))
coeffs={'feature':[],'region':[],'a':[],'err_a':[],'b':[],'err_b':[]}
agg_coeffs={'feature':[],'a':[],'err_a':[],'b':[],'err_b':[]}
for ff,feature in enumerate(scaling_params[1:]):
    
    plt_xy=boeing_data_ww.copy()
    popt,pcov=curve_fit(fit,np.log(plt_xy['resident_pop'].values),np.log(plt_xy[feature].values),p0=(0.0,0.0)) 
    a=popt[0]; 
    err_a=np.sqrt(pcov[0,0])
    #yr_stats=[year,a,err_a]
    b=popt[1]; 
    err_b=np.sqrt(pcov[1,1])
    agg_coeffs['a'].append(a)
    agg_coeffs['err_a'].append(err_a)
    agg_coeffs['b'].append(b)
    agg_coeffs['err_b'].append(err_b)
    agg_coeffs['feature'].append(feature)
    subregions = boeing_data_ww['world_subregion'].drop_duplicates().values
    regions = boeing_data_ww['world_region'].drop_duplicates().values
    list_regions = list(regions)

    for ii,region in enumerate(reversed(list_regions)):
        marker=markers[ii % len(markers)]
        plt_xy=boeing_data_ww.loc[boeing_data_ww['world_region']==region,]
        popt,pcov=curve_fit(fit,np.log(plt_xy['resident_pop'].values),np.log(plt_xy[feature].values),p0=(0.0,0.0)) 
        a=popt[0]; 
        err_a=np.sqrt(pcov[0,0])
        #yr_stats=[year,a,err_a]
        b=popt[1]; 
        err_b=np.sqrt(pcov[1,1])
        coeffs['region'].append(region)
        coeffs['a'].append(a)
        coeffs['err_a'].append(err_a)
        coeffs['b'].append(b)
        coeffs['err_b'].append(err_b)
        coeffs['feature'].append(feature)
        perr = np.sqrt(np.diag(pcov))
        pop = np.array([1,10**9])
        axes[ff].plot(pop,np.exp(b)*pop**a,color=sp((ii-1)/6))
    for ii,region in enumerate(reversed(list_regions)):
        marker=markers[ii % len(markers)]
        plt_xy=boeing_data_ww.loc[boeing_data_ww['world_region']==region,]
        axes[ff].scatter(plt_xy['resident_pop'],plt_xy[feature],color=sp((ii-1)/6),marker=marker,alpha=0.3,label=region)# vary color by world region/country
    axes[ff].set_yscale('log')
    axes[ff].set_xscale('log')
    axes[ff].set_xlim([10**4.5,10**8])
    axes[ff].set_xlabel('Population',**labelfonts)
    if feature == 'length_total':
        axes[ff].set_ylabel('Road Length ($km$)',**labelfonts)
        axes[ff].set_ylim([10**3,10**8.5])
        axes[ff].set_aspect(0.6)
    else:
        if feature == 'area':
            axes[ff].set_ylabel('Urban Center Area ($km^2$)',**labelfonts)
        if feature == 'built_up_area':
            axes[ff].set_ylabel('Builtup Area ($km^2$)',**labelfonts)

        axes[ff].set_ylim([10**-0.1,10**4.1])
        axes[ff].set_aspect(0.8)
    
    axes[ff].tick_params(axis='both', which='major', labelsize=14)
axes[0].legend(loc='lower right', bbox_to_anchor=(1, 0))
plt.savefig('WW_Scaling_region.pdf',transparent=True)
plt.show()
plt.close()



labelfonts = {'fontname':'Arial','fontsize':18}
boeing_data_ww=pd.read_csv('indicators.csv')
country_region = boeing_data_ww[['country_iso','world_subregion']].drop_duplicates()
GDP=pd.read_csv('API_NY.GDP.MKTP.CD_DS2_en_csv_v2_2017804.csv')
Pop=pd.read_csv('API_SP.POP.TOTL_DS2_en_csv_v2_3358390.csv')
Pop.columns = [c+'_pop' if c not in ['Country Name','Country Code','Indicator Name','Indicator Code'] else c for c in Pop.columns]

GDP = pd.merge(left=GDP,right=Pop,on='Country Code')
GDP['country_iso']=GDP['Country Code']
GDP = pd.merge(left=country_region,right=GDP,on='country_iso')

scaling_params = ['resident_pop','area','built_up_area','length_total']
markers=['^','o','s','D','*','v','<','>']#,'1','2','3','4']

fig, ax1 = plt.subplots(1,3,figsize=(17,6))
for ii,feature in enumerate(scaling_params[1:]):
    coeffs=pd.read_csv(feature+'_BoeingWW.csv')
    
    subregion_gdp_ab=[]
    print(coeffs.columns)
    for subregion in GDP['world_subregion'].drop_duplicates():
        region_gdp=0
        num_countries = 0
        for n,row in GDP.loc[GDP['world_subregion']==subregion,].iterrows():
            for year in range(2019,2000,-1):
                latest_gdp = row[str(year)]
                latest_pop = row[str(year)+'_pop']
                if ~np.isnan(latest_gdp) and ~np.isnan(latest_pop):
                    region_gdp+=latest_gdp/latest_pop
                    num_countries += 1
                    break
        region_gdp /= num_countries
        subregion_a=coeffs.loc[coeffs['region']==subregion,'a'].values[0]
        subregion_aerr=coeffs.loc[coeffs['region']==subregion,'err_a'].values[0]
        subregion_b=coeffs.loc[coeffs['region']==subregion,'b'].values[0]
        subregion_berr=coeffs.loc[coeffs['region']==subregion,'err_b'].values[0]
        if subregion!= 'Melanesia':#subregion_a - 2*subregion_aerr > 0:
            subregion_gdp_ab.append([region_gdp,subregion_a,subregion_aerr,subregion_b,subregion_berr])
        else:
            subregion_gdp_ab.append([None]*5)
    subregion_gdp_ab=np.array(subregion_gdp_ab)

    color = 'tab:red'
    theory = {'area':2/3,'built_up_area':2/3,'length_total':5/6}
    coef_label = {'area':'Urban Center Area','built_up_area':'Built-Up Area','length_total':'Road Length'}
    ax1[ii].plot([1,10**7],[theory[feature],theory[feature]],'r--')
    ax1[ii].fill_between([10**10,10**11],[0,0],[1.2,1.2],color='gray',alpha=0.2)
    ax1[ii].set_xlim([0,5.5*10**4])
    ax1[ii].set_ylim([0,1.2])
    
    ax1[ii].set_xlabel('GDP Per Capita',**labelfonts)
    ax1[ii].set_ylabel(coef_label[feature]+' Exponent', color=color,**labelfonts)
    ax1[ii].tick_params(axis='y', labelcolor=color)
    ax1[ii].tick_params(axis='x', rotation=40)
    ax1[ii].tick_params(axis='both', which='major', labelsize=14)
    for jj, region in enumerate(GDP['world_subregion'].drop_duplicates()):
        marker=markers[jj % len(markers)]
        if subregion_gdp_ab[jj,0] is not None:
            ax1[ii].errorbar([subregion_gdp_ab[jj,0]],[subregion_gdp_ab[jj,1]],yerr=[subregion_gdp_ab[jj,2]],color='r',marker='o',markeredgecolor='w',markersize=8,linestyle='',label='Exponent')

    ax2 = ax1[ii].twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel(coef_label[feature]+' Coefficient', color=color,**labelfonts)  # we already handled the x-label with ax1
    for jj, region in enumerate(GDP['world_subregion'].drop_duplicates()):
        marker=markers[jj % len(markers)]
        if subregion_gdp_ab[jj,0] is not None:
            ax2.errorbar([subregion_gdp_ab[jj,0]],[subregion_gdp_ab[jj,3]],yerr=[subregion_gdp_ab[jj,4]],color=color,markersize=8,marker='s',markerfacecolor='none',linestyle='',label='Constant')
    ax2.set_xlim([0,5.5*10**4])
    yvals = subregion_gdp_ab[subregion_gdp_ab[:,0]!=None,3]
    ax2.set_ylim([int(np.min(yvals)-1),int(np.max(yvals)+3)])
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.tick_params(axis='both', which='major', labelsize=14)
    correl_data=subregion_gdp_ab[(subregion_gdp_ab[:,0] != None)]
    ra,pa=spearmanr(list(correl_data[:,0]),list(correl_data[:,1]))
    rc,pc=spearmanr(list(correl_data[:,0]),list(correl_data[:,3]))
plt.tight_layout()
plt.savefig('Const_WorldScaling.pdf')
plt.show()
plt.close()
