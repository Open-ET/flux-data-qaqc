flux-data-qaqc
================

A Python package for correcting/validating the surface energy balance for `FLUXNET climate station data <https://fluxnet.fluxdata.org>`_ (or other climate station data) using techniques that adjust latent energy and sensible heat fluxes. Major functionality includes parsing of climate time series, performing QA/QC routines, and generating visualizations of corrected energy balance components and radiative-turbulent flux ratios at daily and monthly frequencies. The software is currently in early development. 

Documentation
-------------

Online documentation is in development, until then it is recommended to view the `FLUXNET_2015_example <https://github.com/Open-ET/flux-data-qaqc/blob/master/example_data/FLUXNET_2015_example.ipynb>`_ Jupyter notebook. The example includes explanations of multiple functionalities of ``flux-data-qaqc`` and utilizes a provided `FLUXNET 2015 daily climate file <https://raw.githubusercontent.com/Open-ET/flux-data-qaqc/master/examples/FLX_US-AR1_FLUXNET2015_SUBSET_DD_2009-2012_1-3.xlsx>`_.

Installation
------------
Currently clone or download from `GitHub <https://github.com/Open-ET/flux-data-qaqc/edit/master/README.md>`_.  

.. code-block:: bash

    $ git clone https://github.com/Open-ET/flux-data-qaqc.git

Optionally, as opposed to using PIP, you may install dependencies with the provided Conda virtual environment. This is useful to avoid changing your local Python environment. Note, ``flux-data-qaqc`` has been tested for Python 3.7+, although it may work with versions greater than or equal to 3.4.

First make sure you have the ``fluxdataqaqc`` environment file, you can download it `here <https://raw.githubusercontent.com/Open-ET/flux-data-qaqc/master/environment.yml?token=AB3BJKUKL2ELEM7WPLYLXFC45WQOG>`_. Next to install run,

.. code-block:: bash

    $ conda env create -f environment.yml

To activate the environment before using the ``flux-data-qaqc`` package run,

.. code-block:: bash

    $ activate fluxdataqaqc

Run the following to install ``flux-data-qaqc`` in developer mode, soon the package will be uploaded and available on `PYPI <https://pypi.org>`_,

.. code-block:: bash

    $ cd flux-data-qaqc
    $ pip install -e .

Now all package modules and tools should be available in your Python environment PATH and able to be imported. Note if you did not install the Conda virtual environment above, PIP should install dependencies automatically but be sure to be using a version of Python above or equal to 3.4. 
