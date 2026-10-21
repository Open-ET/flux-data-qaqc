flux-data-qaqc
==============

|PyPI| |Downloads| |Tests| |Docs| |JOSS|

``flux-data-qaqc`` is a Python package for energy balance closure correction
and post-processing of eddy covariance flux data. It starts from processed
half-hourly, hourly, or daily station data, such as AmeriFlux and FLUXNET
files, and produces daily and monthly time series of the surface energy
balance components and evapotranspiration (ET) with quality control, gap
filling, and energy balance closure corrections, along with interactive plots.
Applications include hydrologic and atmospheric model validation, evaluation
of satellite-based ET, and irrigation and water use studies.

.. image:: https://raw.githubusercontent.com/Open-ET/flux-data-qaqc/master/docs/source/_static/tutorial/US-Tw3_daily_closure.png
   :alt: Daily energy balance closure at the US-Tw3 AmeriFlux site before and after correction
   :width: 700
   :target: https://flux-data-qaqc.readthedocs.io/en/latest/_static/tutorial/US-Tw3_plots.html

Daily energy balance closure at the US-Tw3 AmeriFlux site before (init) and
after (corr) the Energy Balance Ratio correction, one of the plots made by
``QaQc.plot``. Click the figure to open all of the interactive plots for this
site, or see the
`tutorial <https://flux-data-qaqc.readthedocs.io/en/latest/tutorial.html>`__.

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

Used in
-------

* The eddy covariance ET benchmark dataset for 194 stations in the
  conterminous U.S. (`Volk et al., 2023
  <https://doi.org/10.1016/j.agrformet.2023.109307>`__; `data
  <https://doi.org/10.5281/zenodo.7636781>`__), used to assess the accuracy
  of OpenET (`Volk et al., 2024
  <https://doi.org/10.1038/s44221-023-00181-7>`__).
* Many others

Documentation
-------------

`ReadTheDocs <https://flux-data-qaqc.readthedocs.io/>`_

Installation
------------

Using PIP:

.. code-block:: bash

   pip install fluxdataqaqc

PIP should install the necessary dependencies however it is recommended to use
conda and first install the provided virtual environment. This is useful to
avoid changing your local Python environment. Note, ``flux-data-qaqc``
requires Python 3.9 or newer and is tested with Python 3.9, 3.11, and 3.13.

First make sure you have the ``fluxdataqaqc`` environment file, you can download it `here <https://raw.githubusercontent.com/Open-ET/flux-data-qaqc/master/environment.yml>`_. Next to install run,

.. code-block:: bash

   conda env create -f environment.yml

To activate the environment before using the ``flux-data-qaqc`` package run,

.. code-block:: bash

   conda activate fluxdataqaqc

Now install using PIP:

.. code-block:: bash

   pip install fluxdataqaqc

Now all package modules and tools should be available in your Python environment PATH and able to be imported. Note if you did not install the Conda virtual environment above, PIP should install dependencies automatically but be sure to be using a version of Python above or equal to 3.9. To test that everything has installed correctly by opening a Python interpreter or IDE and run the following:

.. code-block:: python

   import fluxdataqaqc

and 

.. code-block:: python

   from fluxdataqaqc import Data, QaQc, Plot

If everything has been installed correctly you should get no errors. 


How to cite
-----------

Volk et al., (2021). flux-data-qaqc: A Python Package for Energy Balance Closure and Post-Processing of Eddy Flux Data. Journal of Open Source Software, 6(66), 3418, https://doi.org/10.21105/joss.03418


.. |PyPI| image:: https://img.shields.io/pypi/v/fluxdataqaqc.svg
   :target: https://pypi.org/project/fluxdataqaqc/
   :alt: PyPI version

.. |Downloads| image:: https://img.shields.io/pypi/dm/fluxdataqaqc.svg
   :target: https://pypi.org/project/fluxdataqaqc/
   :alt: Downloads per month

.. |Tests| image:: https://github.com/Open-ET/flux-data-qaqc/actions/workflows/fluxdataqaqc_tests.yml/badge.svg
   :target: https://github.com/Open-ET/flux-data-qaqc/actions/workflows/fluxdataqaqc_tests.yml
   :alt: Automated tests

.. |Docs| image:: https://readthedocs.org/projects/flux-data-qaqc/badge/?version=latest
   :target: https://flux-data-qaqc.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation status

.. |JOSS| image:: https://joss.theoj.org/papers/10.21105/joss.03418/status.svg
   :target: https://doi.org/10.21105/joss.03418
   :alt: JOSS paper
