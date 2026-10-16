# -*- coding: utf-8 -*-
"""
Partitioning of net ecosystem exchange (NEE) into gross primary productivity
(GPP) and ecosystem respiration (Reco) with the nighttime method.

Gives GPP estimated the same way at every site from the measured CO2 flux,
e.g. for photosynthetic energy storage (:meth:`.Data.calc_pes`), see
:meth:`.Data.partition_nee`. Low turbulence nights are removed with friction
velocity (u*) thresholds from the moving point test of Papale et al. (2006)
and respiration is modeled with the temperature response of Lloyd and Taylor
(1994) fit to nighttime NEE as in Reichstein et al. (2005). This is similar
to the nighttime partitioning of FLUXNET2015 (ONEFlux) and REddyProc, without
bootstrapping, gap filling or uncertainty estimates, so daytime GPP is only
estimated for records with measured NEE.
"""

import warnings

import numpy as np
import pandas as pd

# Lloyd and Taylor (1994) reference temperature (15 C) and the temperature
# where respiration goes to zero [K]
T_REF = 288.15
T_0 = 227.13


def lloyd_taylor(t_air, r_ref, e_0):
    """
    Ecosystem respiration from air temperature (Lloyd and Taylor, 1994).

    Arguments:
        t_air (array like): air temperature [C].
        r_ref (float or array like): respiration at 15 C, the output has
            the same units.
        e_0 (float or array like): temperature sensitivity [K].

    Returns:
        reco (:obj:`numpy.ndarray`): ecosystem respiration.

    Example:
        Respiration at 15 C is ``r_ref``.

        >>> from fluxdataqaqc.partition import lloyd_taylor
        >>> reco = lloyd_taylor([5, 15, 25], 2, 300)
    """
    t_kelvin = np.asarray(t_air, dtype=float) + 273.15
    return r_ref * np.exp(e_0 * (1 / (T_REF - T_0) - 1 / (t_kelvin - T_0)))


def ustar_thresholds(nee, ustar, t_air, night, n_t_classes=6,
        n_ustar_classes=20, min_records=100, plateau=0.99, max_corr=0.4):
    """
    Friction velocity (u*) threshold of each year with the moving point
    test of Papale et al. (2006), without bootstrapping.

    Nighttime records of each 3 month season (January to March, April to
    June, ...) are split into 6 air temperature classes of equal size, and
    each of these into 20 u* classes of equal size. The threshold of a
    temperature class is the mean u* of the first u* class with mean NEE of
    at least 99% of the mean NEE of the higher u* classes. Temperature
    classes where air temperature and u* are correlated (abs(r) >= 0.4) are
    skipped. The threshold of a season is the median of its temperature
    classes and the threshold of a year the highest of its seasons. Years
    without a threshold get the median of the other years.

    Arguments:
        nee (:obj:`pandas.Series`): sub-daily NEE with a datetime index.
        ustar (:obj:`pandas.Series`): friction velocity [m/s].
        t_air (:obj:`pandas.Series`): air temperature [C].
        night (:obj:`pandas.Series`): True for nighttime records.

    Keyword Arguments:
        n_t_classes (int): default 6. Number of air temperature classes.
        n_ustar_classes (int): default 20. Number of u* classes.
        min_records (int): default 100. Minimum number of records in a
            temperature class.
        plateau (float): default 0.99. Fraction of the mean NEE of the
            higher u* classes that defines the threshold.
        max_corr (float): default 0.4. Temperature classes with an absolute
            correlation of air temperature and u* of this or more are
            skipped.

    Returns:
        thresholds (:obj:`pandas.Series`): u* threshold [m/s] of each year
        (index), null for all years if no year had enough data.
    """
    data = pd.DataFrame({'nee': nee, 'ustar': ustar, 't_air': t_air})
    data = data[night.fillna(False).astype(bool)].dropna()
    years = np.unique(nee.index.year)
    thresholds = pd.Series(np.nan, index=years)

    for year in years:
        in_year = data[data.index.year == year]
        season_thresholds = []
        for quarter in range(1, 5):
            season = in_year[in_year.index.quarter == quarter]
            if len(season) < n_t_classes * min_records:
                continue
            # temperature classes of equal size
            t_class = pd.qcut(
                season.t_air.rank(method='first'), n_t_classes, labels=False
            )
            class_thresholds = []
            for _, rows in season.groupby(t_class):
                # u* that changes with temperature would mix both effects
                corr = np.corrcoef(rows.t_air, rows.ustar)[0, 1]
                if abs(corr) >= max_corr:
                    continue
                threshold = _ustar_class_threshold(
                    rows, n_ustar_classes, plateau
                )
                if threshold is not None:
                    class_thresholds.append(threshold)
            if class_thresholds:
                season_thresholds.append(np.median(class_thresholds))
        if season_thresholds:
            thresholds[year] = max(season_thresholds)

    # years without enough nighttime data use the median of the other years
    thresholds = thresholds.fillna(thresholds.median())
    return thresholds


def _ustar_class_threshold(rows, n_ustar_classes=20, plateau=0.99):
    """
    u* threshold of one temperature class, the mean u* of the first u*
    class whose mean NEE is at least ``plateau`` times the mean NEE of the
    higher u* classes (Papale et al., 2006), None if no class reaches it.
    """
    rows = rows.sort_values('ustar')
    # u* classes of equal size
    u_class = np.arange(len(rows)) * n_ustar_classes // len(rows)
    means = rows.groupby(u_class).mean()
    for i in range(len(means) - 1):
        higher = means.nee.iloc[i + 1:].mean()
        if means.nee.iloc[i] >= plateau * higher:
            return means.ustar.iloc[i]
    return None


def nighttime_partitioning(nee, t_air, night, ustar=None, e0_days=15,
        e0_step=5, rref_days=7, rref_step=4, e0_range=(30, 450),
        min_t_range=5, max_gap=30):
    """
    Partition NEE into GPP and ecosystem respiration with the nighttime
    method of Reichstein et al. (2005).

    Steps:

    1. nighttime NEE with friction velocity below the threshold of its year
       is removed, see :func:`ustar_thresholds`
    2. the temperature sensitivity E0 of each year is the mean of the 3
       short term estimates with the smallest standard errors, from fits of
       :func:`lloyd_taylor` to nighttime NEE in 15 day windows moved by 5
       days, with at least 6 records and a 5 C air temperature range,
       estimates outside 30-450 K are not used
    3. the respiration at 15 C (Rref) is fit with E0 fixed in 7 day windows
       moved by 4 days and linearly interpolated in time, but not across
       gaps of more than 30 days between windows
    4. Reco = Rref f(air temperature, E0), daytime GPP = Reco - NEE and
       GPP = 0 at night (no photosynthesis)

    Years without an E0 estimate use the mean of the other years. Daytime
    GPP is null where NEE or air temperature is missing, like the measured
    (quality flag 0) NEE of FLUXNET2015. Air temperature outside -45 to 60 C
    is treated as missing (the respiration model is not defined at or below
    -46.02 C).

    Arguments:
        nee (:obj:`pandas.Series`): sub-daily NEE with a datetime index,
            positive for a flux to the atmosphere, e.g. [umol m-2 s-1].
        t_air (:obj:`pandas.Series`): air temperature [C].
        night (:obj:`pandas.Series`): True for nighttime records.

    Keyword Arguments:
        ustar (:obj:`pandas.Series` or None): default None. Friction
            velocity [m/s], nighttime data are not u* filtered without it.
        e0_days (int): default 15. Window length for E0 [days].
        e0_step (int): default 5. Step between E0 windows [days].
        rref_days (int): default 7. Window length for Rref [days].
        rref_step (int): default 4. Step between Rref windows [days].
        e0_range (tuple): default (30, 450). Range of valid E0 [K].
        min_t_range (float): default 5. Minimum air temperature range in an
            E0 window [C].
        max_gap (float): default 30. Longest gap between Rref windows that
            is interpolated [days], Reco and GPP are null in longer gaps.

    Returns:
        result (:obj:`pandas.DataFrame`): sub-daily "gpp" and "reco" in the
        units of NEE, and yearly (:obj:`pandas.DataFrame`): u* threshold
        "ustar_threshold" [m/s], "E0" [K] and the number of nighttime
        records used "n_night" of each year.

    Raises:
        ValueError: if E0 could not be estimated in any year.

    Example:
        >>> from fluxdataqaqc.partition import nighttime_partitioning
        >>> night = sw_in < 10
        >>> result, yearly = nighttime_partitioning(
        ...     nee, t_air, night, ustar=ustar
        ... )
        >>> gpp = result.gpp
    """
    night = night.fillna(False).astype(bool)
    # respiration model is not defined near T_0
    t_air = t_air.where((t_air > -45) & (t_air < 60))
    years = np.unique(nee.index.year)
    yearly = pd.DataFrame(
        index=years, columns=['ustar_threshold', 'E0', 'n_night'],
        dtype=float
    )

    # nighttime records for the respiration fits, without low u* records
    use = night & nee.notna() & t_air.notna()
    if ustar is not None:
        thresholds = ustar_thresholds(nee, ustar, t_air, night)
        yearly['ustar_threshold'] = thresholds
        if thresholds.notna().any():
            record_threshold = thresholds.reindex(nee.index.year).values
            use = use & (ustar >= record_threshold)
        else:
            print(
                'WARNING: not enough nighttime data to estimate u* '
                'thresholds, nighttime NEE is not u* filtered.'
            )
    yearly['n_night'] = use.groupby(nee.index.year).sum()
    for year in yearly.index[yearly.n_night == 0]:
        print(
            'WARNING: no nighttime NEE for partitioning in {}, GPP is only '
            'estimated within {} days of other data.'.format(year, max_gap)
        )

    # time in days since the start of the record
    days = (nee.index - nee.index[0].normalize()) / pd.Timedelta(days=1)
    days = pd.Series(days, index=nee.index)

    # temperature sensitivity of each year from short term fits
    for year in years:
        in_year = use & (nee.index.year == year)
        yearly.loc[year, 'E0'] = _year_e0(
            nee[in_year], t_air[in_year], days[in_year], e0_days, e0_step,
            e0_range, min_t_range
        )
    if yearly.E0.isna().all():
        raise ValueError(
            'ERROR: the temperature sensitivity of respiration could not be '
            'estimated from nighttime NEE in any year.'
        )
    yearly['E0'] = yearly.E0.fillna(yearly.E0.mean())

    # temperature response of each record with the E0 of its year
    e_0 = yearly.E0.reindex(nee.index.year).values
    response = pd.Series(lloyd_taylor(t_air, 1, e_0), index=nee.index)

    # reference respiration in short windows, least squares with E0 fixed
    night_days = days[use].values
    night_nee = nee[use].values
    night_response = response[use].values
    centers = []
    r_refs = []
    half = rref_days / 2
    for center in np.arange(half, days.max() + rref_step, rref_step):
        # records are in time order, find the window by searching
        start = np.searchsorted(night_days, center - half)
        end = np.searchsorted(night_days, center + half)
        if end - start < 3:
            continue
        x = night_response[start:end]
        r_ref = np.sum(night_nee[start:end] * x) / np.sum(x ** 2)
        if r_ref > 0:
            centers.append(center)
            r_refs.append(r_ref)
    if not centers:
        raise ValueError(
            'ERROR: reference respiration could not be estimated from '
            'nighttime NEE.'
        )
    r_ref = pd.Series(np.interp(days, centers, r_refs), index=nee.index)
    # not interpolated across long gaps between windows
    centers = np.array(centers)
    i = np.clip(np.searchsorted(centers, days.values), 1, len(centers) - 1)
    nearest = np.minimum(
        np.abs(days.values - centers[i - 1]), np.abs(days.values - centers[i])
    )
    r_ref[nearest > max_gap / 2] = np.nan

    reco = r_ref * response
    # no photosynthesis at night
    gpp = (reco - nee).where(~night)
    gpp[night & reco.notna()] = 0
    result = pd.DataFrame({'gpp': gpp, 'reco': reco})
    return result, yearly


def _year_e0(nee, t_air, days, e0_days=15, e0_step=5, e0_range=(30, 450),
        min_t_range=5, n_best=3):
    """
    Temperature sensitivity E0 [K] of one year, the mean of the ``n_best``
    short term estimates with the smallest standard errors (Reichstein et
    al., 2005), null if there is no valid estimate.
    """
    # scipy is installed with scikit-learn
    from scipy.optimize import OptimizeWarning, curve_fit

    if len(nee) == 0:
        return np.nan
    days = days.values
    estimates = []
    half = e0_days / 2
    for center in np.arange(days.min(), days.max() + e0_step, e0_step):
        # records are in time order, find the window by searching
        start = np.searchsorted(days, center - half)
        end = np.searchsorted(days, center + half)
        t = t_air.values[start:end]
        y = nee.values[start:end]
        if len(y) < 6 or t.max() - t.min() < min_t_range:
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', OptimizeWarning)
                params, cov = curve_fit(
                    lloyd_taylor, t, y, p0=[max(np.mean(y), 0.1), 200],
                    maxfev=2000
                )
        except (RuntimeError, ValueError):
            continue
        e_0 = params[1]
        se = np.sqrt(cov[1, 1]) if np.isfinite(cov[1, 1]) else np.nan
        if e0_range[0] <= e_0 <= e0_range[1] and np.isfinite(se):
            estimates.append((se, e_0))
    if not estimates:
        return np.nan
    estimates.sort()
    return np.mean([e for _, e in estimates[:n_best]])
