
import pandas as pd
import numpy as np

"""
Template workflow to do qaqc for correcting energy balance
will change to be used with different datatypes and sources

working with daily time series named "df" with the following 
variables names:

Rn = 'Net radiation, W/m2'
G = 'Soil-heat flux, W/m2'
LE = 'Latent-heat flux, W/m2'
H = 'Sensible-heat flux, W/m2'

and "et_df" = dataframe containing daily actual ET, this will be removed 
"""

# for plotting later
df['Energy'] = df['Net radiation, W/m2'] - df['Soil-heat flux, W/m2']
df['Flux']   = df['Latent-heat flux, W/m2'] + df['Sensible-heat flux, W/m2']

# new columns for 'Adjusted" LE, H 
df['Adj. LE'] = np.NaN
df['Adj. H']  = np.NaN

# numpy arrays for calculations and graphs
Rn     = np.array(df['Net radiation, W/m2'])
G      = np.array(df['Soil-heat flux, W/m2'])
LE     = np.array(df['Latent-heat flux, W/m2'])
H      = np.array(df['Sensible-heat flux, W/m2'])
energy = np.array(df['Energy'])
flux   = np.array(df['Flux'])

# length of daily time series in for loop
length_sub = len(Rn)

# daily average bowen ratio
df['Bowen ratio'] = (df['Sensible-heat flux, W/m2'] / df['Latent-heat flux, W/m2'])

# daily average bowen ratio used in following adjustments
B = np.array(df['Bowen ratio'])                       

# need to adjust LE, H, and flux and then throw them into these arrays
LE_adj = np.array(df['Adj. LE'])
H_adj  = np.array(df['Adj. H'])

# compute adjusted turbulent fluxes for when Rn > 0
# correcting LE and H, method may be faster as function and vectorized
for i in range(0,length_sub):
    if Rn[i] > 0:
        LE_adj[i] = (Rn[i] - G[i]) / (1 + B[i])
        H_adj[i]  = (B[i] / (1 + B[i])) * (Rn[i] - G[i])
        
    else:
        LE_adj[i] = LE[i]
        H_adj[i]  = H[i]
 
# define new object to fill with corrected turbulent fluxes and corrected ET         
flux_length = len(Rn)
flux_adj = np.full(flux_length, np.NaN)
ET_mm_adj = np.full(flux_length, np.NaN)

# compute adjusted turbulent fluxes for when Rn > 0 and Bowen ratio < 0.05
# instead of forcing closure when bowen ratio is often <- 0.8 or threshold when 
# Rn < 0. The average between i-1 and i+1 is taken. If closure is forced by 
# partitioning the error equally between LE and H, LE is drastically increased
# during periods of possibly "bad" data, usually while the measured LE is going down

for i in range(0, length_sub):
    if Rn[i] > 0 and B[i] < 0.05:
        LE_adj[i] = ((LE[i - 1]) + (LE[i + 1]))/2
        H_adj[i]  = ((H[i - 1])  + (H[i + 1]))/2
    
        # force closure? not in use at the moment
        #LE_adj[i] = LE[i] + 0.5 * (energy[i] - flux[i])
        #H_adj[i]  = H[i]  + 0.5 * (energy[i] - flux[i])
        
        flux_adj[i] = LE_adj[i] + H_adj[i]
        
        # If adjusted fluxes are less than original fluxes, keep originals
        if LE_adj[i] < LE[i]:
            LE_adj[i] = LE[i]
    
        if H_adj[i] < H[i]:
            H_adj[i] = H[i]
    
        # convert LE (W/m2) to ET (mm/day)
        ET_mm_adj[i] = 86400 * (LE_adj[i]/(2500000 * 1000)) * 1000
    
# ET_adj sum for respective water year
et_mm_sum = np.nansum(ET_mm_adj)
    
# append adjusted fluxes and ET to original ET dataframe
et_mm.insert(1,column='ET_adj (mm)',value=ET_mm_adj)
et_mm.insert(2,column='Rn (W/m2)',value=Rn)
et_mm.insert(3,column='LE (W/m2)',value=LE)
et_mm.insert(4,column='LE_adj (W/m2)',value=LE_adj)
et_mm.insert(5,column='H (W/m2)',value=H)
et_mm.insert(6,column='H_adj (W/m2)',value=H_adj)
et_mm.insert(7,column='G (W/m2)',value=G)
et_mm.insert(8,column='ET_adj_total (mm)',value=et_mm_sum)

### EBC comparison
# daily EBC columns 
et_mm['EBC'] = \
    (et_mm['H (W/m2)'] + et_mm['LE (W/m2)']) / \
    (et_mm['Rn (W/m2)'] - et_mm['G (W/m2)'])

et_mm['EBC_adj'] = \
    (et_mm['H_adj (W/m2)'] + et_mm['LE_adj (W/m2)']) / \
    (et_mm['Rn (W/m2)'] - et_mm['G (W/m2)'])

# replace undefined/infinity with nan
et_mm['EBC'] = et_mm['EBC'].replace([np.inf,-np.inf],np.nan)
et_mm['EBC_adj'] = et_mm['EBC_adj'].replace([np.inf,-np.inf],np.nan)

### getting daily sums of adjusted ETa in the next few lines for print statement

# sum subhourly measured and adjusted ET to daily totals
ET_meas_df = pd.DataFrame(et_mm['ET_measured (mm)'])
ET_meas_array = np.array(ET_meas_df)
ET_out_df = pd.DataFrame(et_mm['ET_adj (mm)'])
ET_out_array = np.array(ET_out_df)

# totals of measured vs. adjusted ETa
et_meas_total = np.nansum(ET_meas_array)
et_adj_total  = np.nansum(ET_out_array)
print('Total Water Year measured ETa at %s is %s mm' %(n,et_meas_total))
print('Total Water Year adjusted ETa at %s is %s mm' %(n,et_adj_total))

# print statement of mean EBC for the period 
# EBC result depends on which method (daily average or daytime daily average B) in line 174
EBC_mean     = np.nanmean(et_mm['EBC'])
EBC_mean_adj = np.nanmean(et_mm['EBC_adj'])
print('Mean daily EBC at %s is %s ' %(n,EBC_mean))
print('Mean daily adjusted EBC at %s is %s ' %(n,EBC_mean_adj))

#### scatter plot
### Turbulent flux (LE+H) vs. available energy (Rn-G)

# specs of graph axes and tick marks
[x_min,x_max,x_int] = [-100,301,100]
[y_min,y_max,y_int] = [-100,301,100]
x_ticks = np.arange(x_min,x_max,x_int)    
y_ticks = np.arange(y_min,y_max,y_int)

# 1:1 line
x = np.array([-80,280])
y = np.array([-80,280])
    
          
plt.figure(figsize=(5,3.5),dpi=150)
plt.clf()
plt.plot(energy,flux,marker='x',color='blue',ms=4,linestyle='',label='Raw Data')
plt.plot(energy,flux_adj,marker='o',color='green',fillstyle='none',ms=4.5,linestyle='',label='BR Corrected')
plt.plot(x,y,ls='-',color='black',linewidth=0.5)
plt.legend(loc=0,fontsize=6)
plt.title('Daily Energy Balance Correction for %s' %n)
plt.xlabel('Rn - G (W/m2)',fontsize=7,labelpad=8)
plt.ylabel('LE + H (W/m2)',fontsize=7,labelpad=8)
plt.xticks(x_ticks, fontsize=6)
plt.xlim([x_min, x_max])
plt.yticks(y_ticks, fontsize=6)
plt.ylim([y_min, y_max])
plt.tight_layout()
plt.grid()
plt.show()
file_plot_name = 'EBC_'+n+'.png'
plot_path = os.path.join(work_dir, 'figs', file_plot_name)
plt.savefig(plot_path)

