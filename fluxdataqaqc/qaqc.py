from data import Data
import numpy as np
import sys

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

# Initial setup
# Check if user has passed in a config file, or else just grab the default.
print("\nSystem: starting flux qaqc script.")
if len(sys.argv) == 2:
    config_path = sys.argv[1]
else:
    config_path = 'fluxnet_config.ini'

# Read in climate file and import vars to data frame
data = Data(config_path)

# get length of data set
data_length = len(data.df.stringdate)

# create energy, flux, and bowen ratio columns within dataframe
data.df['energy'] = data.df.net_rad - data.df.g_flux
data.df['flux'] = data.df.le_flux + data.df.h_flux
data.df['bowen_ratio'] = data.df.h_flux / data.df.le_flux

# numpy arrays of dataframe vars
net_rad = np.array(data.df.net_rad)
g_flux = np.array(data.df.g_flux)
le_flux = np.array(data.df.le_flux)
h_flux = np.array(data.df.h_flux)
bowen = np.array(data.df.bowen_ratio)
# numpy arrays of new vars
le_flux_adj = np.full(data_length, np.NaN)
h_flux_adj = np.full(data_length, np.NaN)
flux_adj = np.full(data_length, np.NaN)

# compute adjusted turbulent fluxes for when Rn > 0
# correcting LE and H, method may be faster as function and vectorized
for i in range(0, data_length):
    if net_rad[i] > 0:
        le_flux_adj[i] = (net_rad[i] - g_flux[i]) / (1 + bowen[i])
        h_flux_adj[i] = (bowen[i] / (1 + bowen[i])) * (net_rad[i] - g_flux[i])
        
    else:
        le_flux_adj[i] = le_flux[i]
        h_flux_adj[i] = h_flux[i]

# compute adjusted turbulent fluxes for when Rn > 0 and Bowen ratio < 0.05
# instead of forcing closure when bowen ratio is often <- 0.8 or threshold when 
# Rn < 0. The average between i-1 and i+1 is taken. If closure is forced by 
# partitioning the error equally between LE and H, LE is drastically increased
# during periods of possibly "bad" data, usually while the measured LE is going down

for i in range(0, data_length):
    if net_rad[i] > 0 and bowen[i] < 0.05:
        le_flux_adj[i] = ((le_flux[i - 1]) + (le_flux[i + 1]))/2
        h_flux_adj[i] = ((h_flux[i - 1]) + (h_flux[i + 1]))/2

    flux_adj[i] = le_flux_adj[i] + h_flux_adj[i]
        
    # If adjusted fluxes are less than original fluxes, keep originals
    if le_flux_adj[i] < le_flux[i]:
        le_flux_adj[i] = le_flux[i]
    
    if h_flux_adj[i] < h_flux[i]:
        h_flux_adj[i] = h_flux[i]


# add le_flux_adj, h_flux_adj, and flux_adj to dataframe
# TODO: flux_adj in blake's code is placed to not reflect final h_flux and le_flux adj values, confirm with him
data.df['le_flux_adj'] = le_flux_adj
data.df['h_flux_adj'] = h_flux_adj
data.df['flux_adj'] = flux_adj

# add ET/EBC columns to dataframe using le_flux and h_flux of various sources
#  _reg uses uncorrected le_flux and h_flux
#  _adj uses ajusted le_flux and h_flux from our bowen ratio corrections
#  _corr uses corrected le_flux and h_flux as found in data file (if provided)
data.df['et_reg'] = 86400 * (data.df.le_flux/(2500000 * 1000)) * 1000
data.df['et_adj'] = 86400 * (data.df.le_flux_adj/(2500000 * 1000)) * 1000
data.df['et_corr'] = 86400 * (data.df.le_flux_corr /(2500000 * 1000)) * 1000

data.df['ebc_reg'] = (data.df.h_flux + data.df.le_flux) / (data.df.net_rad - data.df.g_flux)
data.df['ebc_adj'] = (data.df.h_flux_adj + data.df.le_flux_adj) / (data.df.net_rad - data.df.g_flux)
data.df['ebc_corr'] = (data.df.h_flux_corr + data.df.le_flux_corr) / (data.df.net_rad - data.df.g_flux)

# replace undefined/infinity with nans in all EBC columns
data.df.ebc_reg = data.df.ebc_reg.replace([np.inf, -np.inf], np.nan)
data.df.ebc_adj = data.df.ebc_adj.replace([np.inf, -np.inf], np.nan)
data.df.ebc_corr = data.df.ebc_corr.replace([np.inf, -np.inf], np.nan)

# TODO: stopping here for now, all EBC code has been tested and seems to work but has not yet been plotted
# also eventually remove numpy arrays from memory