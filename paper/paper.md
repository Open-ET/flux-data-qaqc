---
title: 'flux-data-qaqc: A Python Package for Energy Balance Closure and Post-Processing of Eddy Flux Data'
tags:
  - Python
  - eddy-covariance
  - energy balance closure
  - post-processing
authors:
  - name: John M. Volk^[corresponding author]
    orcid: 0000-0001-9994-1545
    affiliation: 1
  - name: Justin L. Huntington
    affiliation: 1
  - name: Richard Allen
    affiliation: 2
  - name: Forrest Melton
    affiliation: 3
affiliations:
 - name: Desert Research Institue
   index: 1
 - name: University of Nebraska Lincoln
   index: 2
 - name: National Aeronautics and Space Administration
   index: 3
date: 27 April 2021
bibliography: paper.bib
---

# Introduction

Eddy Covariance (EC) systems provide continuous direct measurements of turbulent and radiative heat and energy fluxes between that atmosphere and land surface [@baldocci1998]. Evapotranspiration (ET) estimates from EC systems are highly valuable for many scientific disciplines, policy makers, and operational applications such as irrigation management. Although the EC technique is considered accurate under certain atmospheric conditions, inherit theoretical and practical limitations of the technique are magnified by data uncertainty from post-processing methods. Error can result from the well observed energy balance closure problem resulting from inherit limitations in the technique but also from data processing error, instrumentation error, and other sources [@foken2008; @stoy2013]. Several approaches of post-processing and corrections of flux ET estimates exist, notably those implemented by major flux measurement networks and datasets like FLUXNET2015 [@pasterello2020]. However, post-processing routines often differ by individual site teams, data networks, and end users, and in some cases post-processing decisions are difficult to reproduce. 

# Statement of need

The ``flux-data-qaqc`` open-source Python package provides a standard and reproducible framework for daily energy balance closure corrections and common post-processing steps of EC and meteorological data. Ultimately, final ET estimates from EC data involves multiple human-driven processing decisions and this presents an uncertainty issue when comparing data that is processed using different methods. ``flux-data-qaqc`` does not address the implementation of EC systems or the initial processing of high frequency EC data of which several software exist, instead, it provides a reproducible framework for common post-processing methods and energy balance closure assessment [@wilson2002]. It is focused towards processing daily and monthly ET from input data is around half-hourly or hourly frequency such as AmeriFlux data. Below is a list of post-processing tools that ``flux-data-qaqc`` provides:

* ability to read a variety of tabular data formats using a configuration file
* tools for batch processing via writing configuration files programmatically
* automatic unit conversions
* filtering of poor quality data (flag-based and numerical)
* gap-filling of hourly or higher frequency EC data with day and night time treated separately 
* temporal aggregation of high frequency data to daily and monthly
* interactive visual inspection tools at pre-processed data frequency and for daily and monthly post-processed data, for both EC and meteorological data
* management and archiving of pre- and post-processed data and metadata, including standardized naming
* energy balance closure correction routines:  an Energy Balance Ratio method, the Bowen Ratio method, and a multiple linear regression method
* other atmospheric calculations, including the hourly and daily the American Society of Civil Engineers standardized reference ET equation [@allen2005]
* ability to download gridded weather data and reference ET that corresponds with the EC tower [@abatzoglou2013]
* daily gap-filling of ET using gridded reference ET and site-based fraction of reference ET 

Online documentation goes into detail of the applications, including a detailed tutorial with usage examples. Easy to follow installation instructions via the Python Package Index and a virtual environment for automatic third party software dependency management are included. 

# Design and features

``flux-data-qaqc`` follows a simple hierarchical structure with a few low-level classes that are then inherited in a linear fashion by two classes that are more specific to EC post-processing (\autoref{fig:fig1}).  The software is also modular to facilitate incorporation of additional classes and methods. It would be straightforward to add additional QAQC or closure routines to the package due to the standardized internal naming schemes and data formats shared by the core classes and their attributes. 

![Diagram showing the relationship between software components, arrows signify the data flow pipeline and class inheritance, and green highlighted boxes indicate software whereas non-highlighted indicate input or output data products of the software.\label{fig:fig1}](figure1.pdf)

Input to ``flux-data-qaqc`` is tabular timeseries EC/meterological data and site metadata which are parsed using a configuration file that can be written by the user or using a tool provided with the package; full example files are provided. The Data class parses the data and has tools for unit conversion, quality-based filtering, visualization, and some atmospheric calculations such as dew point temperature and potential solar radiation. From here, the QaQc class is used for gap-filling, daily and monthly aggregation, energy balance closure correction, writing data using a standard format, and creating interactive diagnostic plots. The Data and QaQc classes leverage the Pandas Python library, giving users access to pre- and post-processed data in a Pandas dataframe which facilitates a myriad of other Python-based scientific analyses via Numpy, SciPy and many other packages. The Plot class contains routines for creating standardized diagnostic plots of pre- and post-processed data, it also provides robust API for generating interactive (via the Bokeh package) scatter and time series plots of arbitrary time series data. The Plot and Convert class along with the Util module are used or inherited throughout the package; these low-level tools are modular enough to easily be used in other software or custom data pipelines. 

# Comparison to other software

``flux-data-qaqc`` is most similar to the ONEFlux pipeline in that it includes the same energy balance closure correction with slight modifications [@pasterello2020]. However, there are many differences, for one ONEFlux is a command line application as opposed to an object oriented framework. ONEFlux has many other uses outside of the scope ``flux-data-qaqc`` which is focused on water vapor flux. Similarly, several tools in ``flux-data-qaqc`` are not included in ONEFlux for example: reading of data in different formats via a configuration file, differing hourly gap-filling methods, a daily gap-filling method using gridded reference ET, visualization tools, atmospheric calculations such as standardized reference ET, monthly aggregation, alternative energy balance closure correction routines, etc.

Other EC data processing and QAQC software can be classified into one of more of these groups: 1) processing of high frequency EC data; and 2) flux partitioning, e.g. partitioning water vapor flux into transpiration and direct evaporation components; 3) friction velocity (ustar) filtering of CO2; 4) post-processed gap-filling and filtering; and 5) flux footprint generation. ``flux-data-qaqc`` is dedicated to tasks in group 4 and has methods and scope unique to other open-source post-processing software known to the authors. Besides ONEFlux, here is a list of all other known and similar EC post-processing software, their purposes and some major differences:
        * MDI Meteo gap-filling tool. Online application that performs gap-filling and uses a different methodology.
        * GaFir. R package that performs gap-filling and uses a different methodology.
        * ReddyProc [@wutzler2018]. R package that performs ustar filtering, gap-filling, and flux partitioning. Differences include: CO2 focused, no closure corrections, different gap-filling methods for half-hourly/hourly, no daily gap-filling, no monthly aggregation or plots, strict input formatting, missing ancillary atmospheric calculations. 
        * openeddy. R package that performs data handling, summarizing, and plotting. Differences include: no closure corrections, different QAQC procedures, no gap-filling, strict input formatting, missing ancillary atmospheric calculations.
* PyFluxPro. Python GUI application that performs standardized post-processing and QAQC routines, gap-filling, ustar filtering, flux partitioning, and plotting. Differences include: different gap-filling methods, no closure corrections, no daily gap-filling, missing ancillary atmospheric calculations.
* hesseflux [@cuntz2020]. Python package similar to ReddyProc with similar differences, namely lack of energy balance closure corrections.

# Research enabled by ``flux-data-qaqc``

This package was designed for and is being used to generate a benchmark ET dataset for ground validation and intercomparison with remote sensing ET models that are part of OpenET (Melton et al. 2020). OpenET is a large collaborative effort to provide satellite based ET estimates at the field scale across the western U.S.A. which will greatly improve water management practices among other uses. This package is also being used for processing daily ET from GrapeX vineyard eddy flux sites to evaluate spectral-based ET approaches (Kustas et al, 2018).

# Acknowledgements

Development of ``flux-data-qaqc`` for OpenET is supported by the S.D. Bechtel, Jr. Foundation; the Gordon and Betty Moore Foundation; the Walton Family Foundation; the Windward Fund; the Water Funder Initiative; the North, Central, and South Delta Water Agencies; the NASA Applied Science Program and the NASA Western Water Applications Office; the USGS Landsat Science Team; the California State University Agricultural Research Institute; and the Idaho Agricultural Experiment Station and Nebraska Agricultural Experiment Station.


