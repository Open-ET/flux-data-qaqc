==========================================================
flux-data-qaqc - Tools for Energy Balance Closure Analysis
==========================================================



|Docs| |Downloads per month| |PyPI version| 

-----------

`View on GitHub <https://github.com/Open-ET/flux-data-qaqc>`_. 

``flux-data-qaqc`` provides a framework to create reproducible workflows for validation and analysis of eddy covariance data. The package is intended for those who need to post-process flux data, particularly for generating daily and monthly evapotransipration (ET) timeseries estimates with energy balance closure corrections applied. Applications where this software may be useful include analysis involving eddy covariance flux tower data, hydrologic or atmospehric model validation, and irrigation and water consumption studies. 

Key functionalities and tools include:

* data validation with methods for quality-based filtering
* time series data tools, e.g. temporal aggregation and resampling
* management of site metadata, data provenance, and file structure
* energy balance closure algorithms and other meterological calculations
* downloading and management of `gridMET <http://www.climatologylab.org/gridmet.html>`__ meterological data
* customizable and interactive data visualizations
* batch processing and unit conversions

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   Configuration Options <advanced_config_options>
   Tutorial <tutorial>
   Closure Algorithms <closure_explanation>
   api
   Tests <software_tests>
   contributors

.. toctree::
   :maxdepth: 1
   
   changelog

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |Docs| image:: https://readthedocs.org/projects/flux-data-qaqc/badge/?version=latest
   :target: https://flux-data-qaqc.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. |Downloads per month| image:: https://img.shields.io/pypi/dm/fluxdataqaqc.svg
   :target: https://pypi.python.org/pypi/fluxdataqaqc/

.. |PyPI version| image:: https://img.shields.io/pypi/v/fluxdataqaqc.svg
   :target: https://pypi.python.org/pypi/fluxdataqaqc/
