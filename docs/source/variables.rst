Variables and Units
===================

Calculated variables
--------------------

Names, descriptions, and units of the variables that ``flux-data-qaqc``
calculates or renames using its standardized naming scheme. The names in the
"Variable" column are the default column names in the daily and monthly
output CSV files and are shown in the interactive plots when hovering with a
cursor. To write output files with the names from your input file instead,
for variables that were read in and not modified, use
``use_input_names=True`` with :meth:`.QaQc.write`.

Not all of these variables exist in every case, which ones are in the output
depends on the input data provided in the configuration file and on which
calculations are used. For example, ASCE_ETo will only exist if the hourly or
daily methods for computing ASCE standardized reference ET are used and the
required input variables exist. Variables created by the energy balance
closure corrections are explained in :ref:`Closure Methodologies` and how
input names are mapped to these names is explained in :ref:`Variable names
and units`.


 ================= ================================================================================ =======================
 Variable           Description                                                                     Unit
 ================= ================================================================================ =======================
 ASCE_ETo           short/grass ASCE standardized reference ET                                      mm time⁻¹
 ASCE_ETr           tall/alfalfa ASCE standardized reference ET                                     mm time⁻¹
 blh                boundary layer height                                                           m
 br                 bowen ratio                                                                     —
 co2                CO2 mole fraction                                                               μmol mol⁻¹
 ebc_cf             energy balance closure correction factor (inverse of ebr_est)                   —
 ebc_cf_method      method used for ebc_cf, 1: median of +/- 7 days, 2: mean of +/- 5 days,         —
                    3: +/- 5 days in previous and next years (1-3 follow FLUXNET2015,
                    `Pastorello et al., 2020`_), 4: 5 day climatology (not in FLUXNET)
 ebr                input energy balance ratio                                                      —
 ebr_5day_clim      5 day climatology of the filtered energy balance ratio (method 4)               —
 ebr_corr           energy balance ratio after correction, (H_corr + LE_corr) / (Rn - G)            —
 ebr_est            energy balance ratio estimated by the 'ebr' correction (inverse of ebc_cf)      —
 energy             input Rn - G                                                                    W m⁻²
 es                 saturation vapor pressure                                                       kPa
 ET                 ET calculated from input LE and average air temperature                         mm time⁻¹
 ET_corr            ET calculated from LE_corr                                                      mm time⁻¹
 ET_fill            gridMET_ETr * ETrF_filtered or gridMET_ETo * ETF_filtered (fills gaps in ET)    mm time⁻¹
 ET_gap             True on gap days in ET_corr, False otherwise (n gap-days in monthly files)      —
 ET_user_corr       corrected ET, user-provided                                                     mm time⁻¹
 EToF               fraction of reference ET for ET_corr, i.e. ET_corr / gridMET_ETo                —
 EToF_filtered      filtered and gap-filled EToF                                                    —
 ETrF               fraction of reference ET for ET_corr, i.e. ET_corr / gridMET_ETr                —
 ETrF_filtered      filtered and gap-filled ETrF                                                    —
 fc                 turbulent CO2 flux                                                              μmol CO₂ m⁻² s⁻¹
 flux               input LE + H                                                                    W m⁻²
 flux_corr          LE_corr + H_corr                                                                W m⁻²
 G                  average or single soil heat flux                                                W m⁻²
 G_[1,2,3,…]        soil heat flux at sensor                                                        W m⁻²
 G_subday_gaps      number of gaps in initial G per day                                             —
 gpp                gross primary productivity                                                      μmol CO₂ m⁻² s⁻¹
 gridMET_ETo        gridMET short/grass reference ET (nearest cell)                                 mm time⁻¹
 gridMET_ETr        gridMET tall/alfalfa reference ET (nearest cell)                                mm time⁻¹
 gridMET_prcp       gridMET precipitation (nearest cell)                                            mm time⁻¹
 gridMET_[other]    other optional gridMET variables (nearest cell)                                 NA
 H_corr             corrected sensible heat                                                         W m⁻²
 H_subday_gaps      number of gaps in initial H per day                                             —
 H_user_corr        corrected sensible heat, user-provided                                          W m⁻²
 LE_corr            corrected latent energy                                                         W m⁻²
 LE_subday_gaps     number of gaps in initial LE per day                                            —
 LE_user_corr       corrected latent energy, user-provided                                          W m⁻²
 lw_in              incoming longwave radiation                                                     W m⁻²
 lw_out             outgoing longwave radiation                                                     W m⁻²
 mo_length          Monin-Obukhov length                                                            m
 nee                net ecosystem exchange                                                          μmol CO₂ m⁻² s⁻¹
 pes                photosynthetic energy fixed over the timestep                                   J m⁻²
 pes_flux           photosynthetic energy storage flux-equivalent term                              W m⁻²
 ppfd_in            incoming photosynthetic photon flux density                                     μmol photons m⁻² s⁻¹
 ppt                precipitation                                                                   mm time⁻¹
 reco               ecosystem respiration                                                           μmol CO₂ m⁻² s⁻¹
 rh                 relative humidity                                                               %
 Rn                 net radiation                                                                   W m⁻²
 Rn_subday_gaps     number of gaps in initial Rn per day                                            —
 rso                clear sky radiation (ASCE formulation)                                          W m⁻²
 sc                 CO2 storage flux                                                                μmol CO₂ m⁻² s⁻¹
 sigmav             standard deviation of lateral velocity fluctuations                             m s⁻¹
 sw_in              incoming shortwave radiation                                                    W m⁻²
 sw_out             outgoing shortwave radiation                                                    W m⁻²
 sw_pot             potential shortwave radiation, user-provided                                    W m⁻²
 t_avg              average temperature                                                             C
 t_dew              dew point temperature                                                           C
 t_max              maximum temperature                                                             C
 t_min              minimum temperature                                                             C
 theta              average or single soil moisture                                                 user defined
 theta_[1,2,3,…]    soil moisture at sensor                                                         user defined
 ustar              friction velocity                                                               m s⁻¹
 vp                 actual vapor pressure                                                           kPa
 vpd                vapor pressure deficit                                                          kPa
 wd                 wind direction                                                                  user defined
 ws                 wind speed                                                                      m s⁻¹
 zeta               Monin-Obukhov stability parameter (z/L)                                         —
 ================= ================================================================================ =======================

.. _Pastorello et al., 2020: https://doi.org/10.1038/s41597-020-0534-3


Input units
-----------

Units of input variables are given in the DATA section of the config file,
using the name of the variable's column entry with ``_units`` instead of
``_col``, e.g. ``latent_heat_flux_units = w/m2``. Upper or lower case can be
used. The variables below are used in calculations, so their units must be
given in the config and they are converted to the units in the "Converted
to" column when a :obj:`.QaQc` object is created (and by :obj:`.Data`
methods that use them, e.g. reference ET). Units of other input variables,
e.g. relative humidity, wind direction, and soil moisture, are only used for
plot labels and are not converted.


 ========== ================================================ ====================================== ====================
 Variable   Config units name                                Accepted input units                   Converted to
 ========== ================================================ ====================================== ====================
 LE         ``latent_heat_flux_units``                       ``w/m2``, ``mj/m2``                    ``w/m2``
 H          ``sensible_heat_flux_units``                     ``w/m2``, ``mj/m2``                    ``w/m2``
 Rn         ``net_radiation_units``                          ``w/m2``, ``mj/m2``                    ``w/m2``
 G          ``ground_flux_units``                            ``w/m2``, ``mj/m2``                    ``w/m2``
 lw_in      ``longwave_in_units``                            ``w/m2``, ``mj/m2``                    ``w/m2``
 lw_out     ``longwave_out_units``                           ``w/m2``, ``mj/m2``                    ``w/m2``
 sw_in      ``shortwave_in_units``                           ``w/m2``                               ``w/m2``
 sw_out     ``shortwave_out_units``                          ``w/m2``, ``mj/m2``                    ``w/m2``
 ppt        ``precip_units``                                 ``mm``, ``in``, ``m``                  ``mm``
 vp         ``vap_press_units``                              ``kpa``, ``hpa``, ``pa``               ``kpa``
 vpd        ``vap_press_def_units``                          ``kpa``, ``hpa``, ``pa``               ``kpa``
 t_avg      ``avg_temp_units``                               ``c``, ``f``, ``k``                    ``c``
 ws         ``wind_spd_units``                               ``m/s``, ``mph``                       ``m/s``
 gpp        ``gross_primary_productivity_units``             ``umolco2/m2/s``                       ``umolco2/m2/s``
 reco       ``ecosystem_respiration_units``                  ``umolco2/m2/s``                       ``umolco2/m2/s``
 nee        ``net_ecosystem_exchange_units``                 ``umolco2/m2/s``                       ``umolco2/m2/s``
 fc         ``co2_turbulent_flux_units``                     ``umolco2/m2/s``                       ``umolco2/m2/s``
 sc         ``co2_storage_flux_units``                       ``umolco2/m2/s``                       ``umolco2/m2/s``
 ppfd_in    ``photosynthetic_photon_flux_density_in_units``  ``umolphoton/m2/s``                    ``umolphoton/m2/s``
 co2        ``co2_mole_fraction_units``                      ``umol/mol``, ``ppm``                  ``umol/mol``
 ustar      ``friction_velocity_units``                      ``m/s``                                ``m/s``
 mo_length  ``monin_obukhov_length_units``                   ``m``                                  ``m``
 zeta       ``monin_obukhov_stability_parameter_units``      ``dimensionless``, ``nondimensional``  ``dimensionless``
 sigmav     ``lateral_velocity_fluctuation_std_dev_units``   ``m/s``                                ``m/s``
 blh        ``boundary_layer_height_units``                  ``m``                                  ``m``
 ========== ================================================ ====================================== ====================


.. note::
   ``mj/m2`` means MJ m⁻² per day and is only accepted for daily input data.
   For hourly or shorter input data an error is raised, convert those
   variables to ``w/m2`` first.

The same lists are in the :attr:`.QaQc.allowable_units` and
:attr:`.QaQc.required_units` attributes. The list of allowable units is a
work in progress, if your input units are not available consider raising an
issue on `GitHub <https://github.com/Open-ET/flux-data-qaqc/issues>`__ or
providing the conversion directly with a pull request. Automatic unit
conversions are handled by the :meth:`.util.Convert.convert` class method.
