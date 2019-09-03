==========================================================
flux-data-qaqc - Tools for Energy Balance Closure Analysis
==========================================================

.. image:: https://readthedocs.org/projects/flux-data-qaqc/badge/?version=latest
   :target: https://flux-data-qaqc.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

``flux-data-qaqc`` provides a framework to create reproducible workflows for the analysis of eddy covariance time series data. In particular for performing energy balance closure analysis and correction routines which adjust turbulent fluxes.

Notable tools include:

* data validation with methods for quality-based filtering
* time series data tools, e.g. temporal aggregation and resampling
* management of site metadata, data provenance, and file structure
* energy balance closure algorithms and other meterological calculations
* downloading and management of `gridMET <http://www.climatologylab.org/gridmet.html>`__ meterological data
* customizable and interactive data visualizations
* batch processing
* unit conversions



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   Configuration Options <advanced_config_options>
   Closure Algorithms <closure_explanation>
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
