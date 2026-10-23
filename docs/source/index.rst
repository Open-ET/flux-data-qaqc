==========================================================
flux-data-qaqc - Tools for Energy Balance Closure Analysis
==========================================================



|Docs| |Downloads per month| |PyPI version| 

-----------

`View on GitHub <https://github.com/Open-ET/flux-data-qaqc>`_. 

``flux-data-qaqc`` provides a framework to create reproducible workflows for validation and analysis of eddy covariance data. The package is intended for those who need to post-process flux data, particularly for generating daily and monthly evapotranspiration (ET) timeseries estimates with energy balance closure corrections applied. Applications where this software may be useful include analysis involving eddy covariance flux tower data, hydrologic or atmospheric model validation, and irrigation and water consumption studies. 

Key functionalities and tools include:

* reading sub-daily or daily station data (e.g. AmeriFlux, FLUXNET, or other
  CSV and Excel files) with a config file that holds the site metadata and
  column names
* quality control with QC flags, averaging of multiple sensors, and a filter
  for days with too many sub-daily gaps
* daily and monthly aggregation, with interpolation of short sub-daily gaps
* energy balance closure corrections (energy balance ratio, Bowen ratio, and
  linear regression), optionally with photosynthetic energy storage removed
  from Rn - G
* download of `gridMET <https://www.climatologylab.org/gridmet.html>`__ data
  and gap filling of daily ET with the gridMET reference ET fraction (ETrF)
* partitioning of NEE into GPP and ecosystem respiration
* ASCE reference ET and other meteorological calculations
* interactive plots, unit conversions, and writing config files for batches
  of sites

.. rubric:: Quick start

.. code-block:: bash

   pip install fluxdataqaqc

.. code-block:: python

   from fluxdataqaqc import Data, QaQc

   # read the config file and the station data
   d = Data('US-Tw3_config.ini')
   # resample to daily and monthly, filter days with data gaps
   q = QaQc(d)
   # energy balance closure correction, ET gaps filled with gridMET
   q.correct_data(meth='ebr')
   # save interactive plots (HTML) and daily and monthly CSV files
   q.plot()
   q.write()

Outputs are saved in an "output" folder next to the config file.
The example config file and data are in
`examples/Basic_usage <https://github.com/Open-ET/flux-data-qaqc/tree/master/examples/Basic_usage>`__,
see :ref:`Installation` for other ways to install and the :ref:`Tutorial`
for each step and its options.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   Configuration Options <advanced_config_options>
   Tutorial <tutorial>
   Closure Methodologies <closure_explanation>
   Variables and Units <variables>
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
