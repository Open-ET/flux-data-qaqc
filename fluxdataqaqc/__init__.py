"""
Tools for validating `FLUXNET <https://fluxnet.fluxdata.org>`_ (or other) eddy covariance flux station data. Includes tools for adjusting near surface latent energy and sensible heat fluxes to improve closure of the energy balance using multiple methodologies. Other functionalities include data processing, cleaning data using quality control metrics, and visualization.
"""

__name__ = 'flux-data-qaqc'
__author__ = 'John Volk and Christian Dunkerly'
__version__ = '0.0.1'


from fluxdataqaqc.data import Data
from fluxdataqaqc.qaqc import QaQc
from fluxdataqaqc.plot import Plot
