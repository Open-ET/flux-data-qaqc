"""
Tools for performing energy balance correction and QA/QC of eddy covariance flux station data. For example, routines for adjusting near surface latent energy and sensible heat fluxes to improve closure of the energy balance. Other functionalities include data processing e.g. temporal resampling, interpolation, and aggregation, unit conversions, cleaning data using quality control metrics, linear regression, meterological calcuations, and interactive visualizations.
"""

__name__ = 'fluxdataqaqc'
__author__ = 'John Volk'
__version__ = '0.1.6'


from fluxdataqaqc.data import Data
from fluxdataqaqc.qaqc import QaQc
from fluxdataqaqc.plot import Plot
from fluxdataqaqc.util import Convert
from fluxdataqaqc import util 
