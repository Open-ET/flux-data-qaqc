.. image:: https://readthedocs.org/projects/flux-data-qaqc/badge/?version=latest
   :target: https://flux-data-qaqc.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status


flux-data-qaqc
================

``flux-data-qaqc`` provides a framework to create reproducible workflows for the validation and analysis of eddy covariance time series data.

Notable uses:

* data validation with methods for quality-based filtering
* time series tools, e.g. temporal aggregation, interpolation, and resampling
* energy balance closure algorithms and other meterological calculations
* data provenance, e.g. from metadata management and file structure
* downloading and management of `gridMET <http://www.climatologylab.org/gridmet.html>`__ meterological data
* customizable and interactive visualizations
* built-in unit conversions
* facilitating batch processing 

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
avoid changing your local Python environment. Note, ``flux-data-qaqc`` has been
tested for Python 3.7+, although it may work with versions greater than or
equal to 3.4.

First make sure you have the ``fluxdataqaqc`` environment file, you can download it `here <https://raw.githubusercontent.com/Open-ET/flux-data-qaqc/master/environment.yml?token=AB3BJKUKL2ELEM7WPLYLXFC45WQOG>`_. Next to install run,

.. code-block:: bash

   conda env create -f environment.yml

To activate the environment before using the ``flux-data-qaqc`` package run,

.. code-block:: bash

   conda activate fluxdataqaqc

Now install using PIP:

.. code-block:: bash

   pip install fluxdataqaqc

Now all package modules and tools should be available in your Python environment PATH and able to be imported. Note if you did not install the Conda virtual environment above, PIP should install dependencies automatically but be sure to be using a version of Python above or equal to 3.4. To test that everything has installed correctly by opening a Python interpretor or IDE and run the following:

.. code-block:: python

   import fluxdataqaqc

and 

.. code-block:: python

   from fluxdataqaqc import Data, QaQc, Plot

If everything has been installed correctly you should get no errors. 

