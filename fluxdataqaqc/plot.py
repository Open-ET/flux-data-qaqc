# -*- coding: utf-8 -*-
import re
import numpy as np
from bokeh import models
from bokeh.palettes import Viridis256
from bokeh.layouts import gridplot
from bokeh.plotting import figure, ColumnDataSource, output_file, show, save
from bokeh.models import HoverTool, Label
from bokeh.models.formatters import DatetimeTickFormatter
from pathlib import Path

class Plot(object):
    """
    Container of plot routines of :mod:`fluxdataqaqc` including static methods 
    that can be used to create  and update line and scatter plots from a 
    :obj:`pandas.DataFrame` instance.
    """

    def __init__(self):
        pass

    @staticmethod
    def line_plot(fig, x, y, source, color, label=None,
            x_axis_type='date', **kwargs):
        if label is None:
            fig.line(x,y, source=source, color=color, **kwargs)
        else:
            label=dict(value=label)
            fig.line(x,y, source=source, color=color,legend=label, **kwargs)
            fig.legend.location = "top_left"
        # add Hover tool with additional tips if more than one line
        if len(fig.hover) == 0:
            Hover = HoverTool(
                    tooltips=[
                        (x,'@{}{}'.format(x,'{%F}')), 
                        (y,'@{}'.format(y))
                    ], formatters={x: 'datetime'}
            )
            fig.add_tools(Hover)
        else:
            fig.hover[0].tooltips.append((y,'@{}'.format(y)))
        # enforce datetime x axis it date
        if x_axis_type == 'date':
            fig.xaxis.formatter = DatetimeTickFormatter(days="%d-%b-%Y")


    @staticmethod
    def scatter_plot(fig, x, y, source, color, label='', 
            lsrl=True, date_name='date', **kwargs):
        """
        
        Note:
            If sending extra keyword arguments (for Bokeh) they will be passed
            to the scatter plot but not to the least squares regression 
            line plot.
        """
        name = '{}_vs_{}'.format(x,y)
        # remove pairs where one or both are nans, for LSRL and min/max-1:1
        xd = source.data[x].astype(float)
        yd = source.data[y].astype(float)
        mask = (~np.isnan(xd) & ~np.isnan(yd))
        xd = xd[mask]
        yd = yd[mask]
        # ordinary least squares linear regression slope through zero
        if lsrl:
            #source.data.update(ColumnDataSource({x: xd, y: yd}).data)
            # plot scatter and linear regression line slope for y=mx
            m = np.linalg.lstsq(xd.reshape(-1,1),yd,rcond=None)[0][0]
            fig.circle(
                x, y, source=source, color=color, line_width=1, fill_alpha=0.2,
                legend='{lab}, slope={s:.2f}'.format(lab=label, s=m), name=name,
                size=10, **kwargs
            )
            fig.line(xd, m * xd, color=color)
        else:
            fig.circle(
                x, y, source=source, color=color, line_width=1, fill_alpha=0.5,
                legend=dict(value=label), name=name, size=10, **kwargs
            )
        if len(fig.hover) == 0:
            Hover = HoverTool(
                tooltips=[
                    (date_name,'@{}{}'.format(date_name,'{%F}')),
                    (x,'@{}'.format(x)), 
                    (y,'@{}'.format(y))
                ],formatters={date_name: 'datetime'}
            )
            fig.add_tools(Hover)
        else:
            new_tips = [(x,'@{}'.format(x)),(y,'@{}'.format(y))]
            for t in new_tips:
                if not t in fig.hover[0].tooltips:
                    fig.hover[0].tooltips.append(t)
        fig.hover[0].names.append(name)
        fig.legend.location = "top_left"

        return xd.min(), xd.max()

    @staticmethod
    def add_lines(fig, df, plt_vars, colors, x_name, source, 
            labels=None, **kwargs): 
        """add multiple lines to line plot, if none exist return None"""
        ret = None
        n_lines = 0
        for i, v in enumerate(plt_vars):
            if not v in df.columns or df[v].isna().all():
                continue
            else:
                if labels is None:
                    label = None
                else: 
                    label = labels[i]
                n_lines += 1
                Plot.line_plot(
                    fig, x_name, v, source, color=colors[i], label=label, 
                    **kwargs
                )
        if n_lines > 0:
            ret = fig

        return ret

    def _plot(self, QaQc, ncols=1, output_type='save', out_file=None): 
        """ 
        Private routine for aggregated validation plots which are called by
        default using the :meth:`fluxdataqaqc.QaQc.plot` method. Full API
        documentation is found in the :meth:`fluxdataqaqc.QaQc.plot` method.  
        """
        # get daily and monthly time series with internal names, get units
        df = QaQc.df.rename(columns=QaQc.inv_map) 
        monthly_df = QaQc.monthly_df.rename(columns=QaQc.inv_map) 
        variables = QaQc.variables
        units = QaQc.units 
        if out_file is None:
            out_file = Path(QaQc.out_dir)/'{}_plots.html'.format(QaQc.site_id)
        # add else statement to allow making any subdir that does not yet exist
        # if out_file is to a non-existent directory
        else:
            out_dir = Path(out_file).parent
            if not out_dir.is_dir():
                out_dir.mkdir(parents=True, exist_ok=True)

        output_file(out_file)

        # bokeh column sources for tooltips
        daily_source=ColumnDataSource(df)
        monthly_source = ColumnDataSource(monthly_df)
        # for aggregating plots
        daily_figs = []
        monthly_figs = []

        # run through each plot, daily then monthly versions
        def _get_units(plt_vars, units):
            """
            Helper function to figure out units for multivariate plots.
            If none of plt_vars exist return None, if multiple units are found
            print a warning that vars have different units. Returns string if 
            one or more units are found- first found if multiple. 
            """
            ret = [] 
            for v in plt_vars:
                unit = units.get(v, None)
                if unit is not None:
                    ret.append(unit)
            if len(ret) == 0:
                ret = None
            elif len(set(ret)) > 1:
                print(
                    'WARNING: variables: {} are not of the same units'.format(
                        ','.join(plt_vars)
                    )
                )
                ret = ret[0]
            elif len(set(ret)) == 1:
                ret = ret[0]

            return ret


        #### 
        # energy balance time series plots
        #### 
        plt_vars = ['LE', 'H', 'Rn', 'G']
        colors = ['blue', 'red', 'black', 'green']
        title = 'Daily Surface Energy Balance Components'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(x_axis_label=x_label, y_axis_label=y_label, title=title)
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=plt_vars
        )
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Surface Energy Balance Components'
            fig = figure(x_axis_label=x_label, y_axis_label=y_label,title=title)
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=plt_vars
            )
            monthly_figs.append(fig)
        else:
            print(
                'Energy balance components time series grapths missing all '
                'variables'
            )

        #### 
        # incoming shortwave and ASCE potential clear sky time series plots
        #### 
        plt_vars = ['sw_in', 'rso']
        labels = ['station', 'ASCE_pot']
        colors = ['black', 'red']
        title = 'Daily Incoming Shortwave and ASCE Potential Radiation'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(x_axis_label=x_label, y_axis_label=y_label, title=title)
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=labels
        )
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Incoming Shortwave and ASCE Potential Radiation'
            fig = figure(x_axis_label=x_label, y_axis_label=y_label,title=title)
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=labels
            )
            monthly_figs.append(fig)
        else:
            print(
                'Shortwave and potential clear sky shortwave time series '
                'grapths missing all variables'
            )

        #### 
        # radiation time series plots
        #### 
        plt_vars = ['sw_in', 'lw_in', 'sw_out', 'lw_out']
        colors = ['red', 'darkred', 'blue', 'navy']
        title = 'Daily Radiation Components'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(x_axis_label=x_label, y_axis_label=y_label, title=title)
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=plt_vars
        )
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Radiation Components'
            fig = figure(x_axis_label=x_label, y_axis_label=y_label,title=title)
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=plt_vars
            )
            monthly_figs.append(fig)
        else:
            print(
                'Radiation components time series grapths missing all variables'
            )

        #### 
        # average temperature time series plot
        #### 
        plt_vars = ['t_avg']
        colors = ['black']
        title = 'Daily Average Air Temperature'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(x_axis_label=x_label, y_axis_label=y_label, title=title)
        fig = Plot.add_lines(fig, df, plt_vars, colors, x_label, daily_source)
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Average Air Temperature'
            fig = figure(x_axis_label=x_label, y_axis_label=y_label,title=title)
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source
            )
            monthly_figs.append(fig)
        else:
            print(
                'Average air temperature time series grapths missing all '
                'variables'
            )

        #### 
        # vapor pressure time series plots
        #### 
        plt_vars = ['vp', 'vpd']
        colors = ['blue', 'black']
        title = 'Daily Average Vapor Pressure and Deficit'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(x_axis_label=x_label, y_axis_label=y_label, title=title)
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=plt_vars
        )
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Average Vapor Pressure'
            fig = figure(x_axis_label=x_label, y_axis_label=y_label,title=title)
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=plt_vars
            )
            monthly_figs.append(fig)
        else:
            print('Vapor pressure time series grapths missing all variables')

        #### 
        # windpseed time series plot
        #### 
        plt_vars = ['ws']
        colors = ['black']
        title = 'Daily Average Windspeed'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(x_axis_label=x_label, y_axis_label=y_label, title=title)
        fig = Plot.add_lines(fig, df, plt_vars, colors, x_label, daily_source)
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Average Windspeed'
            fig = figure(x_axis_label=x_label, y_axis_label=y_label,title=title)
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source
            )
            monthly_figs.append(fig)
        else:
            print('Windspeed time series grapths missing all variables')

        #### 
        # precipitation time series plots
        #### 
        plt_vars = ['ppt', 'gridMET_prcp_mm']
        labels = ['station', 'gridMET']
        colors = ['black', 'red']
        title = 'Daily Station and gridMET Precipitation'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(x_axis_label=x_label, y_axis_label=y_label, title=title)
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=labels
        )
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Station and gridMET precipitation'
            fig = figure(x_axis_label=x_label, y_axis_label=y_label,title=title)
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=labels
            )
            monthly_figs.append(fig)
        else:
            print('Precipitation time series grapths missing all variables')

        #### 
        # latent energy time series plots
        #### 
        plt_vars = ['LE', 'LE_corr', 'LE_user_corr']
        colors = ['black', 'red', 'darkred']
        title = 'Daily Average Latent Energy Flux'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(x_axis_label=x_label, y_axis_label=y_label, title=title)
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=plt_vars
        )
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Average Vapor Pressure'
            fig = figure(x_axis_label=x_label, y_axis_label=y_label,title=title)
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=plt_vars
            )
            monthly_figs.append(fig)
        else:
            print('Latent energy time series grapths missing all variables')

        #### 
        # ET time series plots
        #### 
        plt_vars = ['et', 'et_corr', 'et_user_corr', 'et_fill_val']
        colors = ['black', 'red', 'darkred', 'green']
        title = 'Daily Evapotranspiration'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(x_axis_label=x_label, y_axis_label=y_label, title=title)
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=plt_vars
        )
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Evapotranspiration'
            fig = figure(x_axis_label=x_label, y_axis_label=y_label,title=title)
            # do not sum et_fill_val (gap filled with etr*kc) instead ndays fill
            plt_vars[3] = 'et_gap'
            labels = plt_vars[0:3] + ['days_filled']
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=labels
            )
            monthly_figs.append(fig)
        else:
            print(
                'Evapotranspiration time series grapths missing all variables'
            )

        #### 
        # Kc time series plots
        #### 
        plt_vars = ['Kc', 'Kc_7day_mean']
        colors = ['black', 'red']
        title = 'Daily Crop Coefficient (Kc)'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(x_axis_label=x_label, y_axis_label=y_label, title=title)
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=plt_vars
        )
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Crop Coefficient (Kc)'
            fig = figure(x_axis_label=x_label, y_axis_label=y_label,title=title)
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=plt_vars
            )
            monthly_figs.append(fig)
        else:
            print('Crop coefficient time series grapths missing all variables')

        #### 
        # energy balance ratio time series plots
        #### 
        plt_vars = ['ebr', 'ebr_corr', 'ebr_user_corr']
        colors = ['black', 'red', 'darkred']
        title = 'Daily Energy Balance Ratio with Long-term Mean'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        # add mean EBR for each time series in legend
        labels = []
        for i, v in enumerate(plt_vars):
            if v in df.columns:
                added_text = ': {}'.format(str(round(df[v].mean(),2)))
                labels.append(plt_vars[i] + added_text)
            else:
                labels.append(None)
        fig = figure(x_axis_label=x_label, y_axis_label=y_label, title=title)
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=labels
        )
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Energy Balance Ratio with Long-term Mean'
            # add mean for monthly EBRs to legend
            labels = []
            for i, v in enumerate(plt_vars):
                if v in df.columns:
                    added_text = ': {}'.format(
                        str(round(monthly_df[v].mean(),2))
                    )
                    labels.append(plt_vars[i] + added_text)
                else:
                    labels.append(None)
            fig = figure(x_axis_label=x_label, y_axis_label=y_label,title=title)
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=labels
            )
            monthly_figs.append(fig)
        else:
            print(
                'Energy balance ratio time series grapths missing all '
                'variables'
            )

        #### 
        # energy balance closure scatter plots
        #### 
        title = 'Daily Energy Balance Closure, Energy versus Flux with Slope '\
            'Through Zero'
        unit = _get_units(['LE', 'H', 'Rn', 'G'], units)
        y_label = 'LE + H ({})'.format(unit)
        x_label = 'Rn - G ({})'.format(unit)
        fig = figure(x_axis_label=x_label, y_axis_label=y_label, title=title)
        y_vars = ['flux', 'flux_corr', 'flux_user_corr']
        colors = ['black', 'red', 'darkred']
        labels = ['init', 'corr', 'user_corr']
        # add plot pairs to plot if they exist, add 1:1
        mins_maxs = []
        n_vars_fnd = 0
        for i, v in enumerate(y_vars):
            if v in df.columns and not df[v].isna().all():
                n_vars_fnd += 1
                min_max = Plot.scatter_plot(
                    fig, 'energy', v, daily_source, colors[i], label=labels[i]
                )
                mins_maxs.append(min_max)
        if n_vars_fnd > 0:
            # add scaled one to one line
            mins_maxs = np.array(mins_maxs)
            one2one_vals = np.arange(min(mins_maxs[:,0]), max(mins_maxs[:,1]),1)
            fig.line(
                one2one_vals, one2one_vals, legend='1:1 line', color='black', 
                line_dash='dashed'
            )
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Energy Balance Closure, Energy versus Flux '\
                'with Slope Through Zero'
            fig = figure(x_axis_label=x_label, y_axis_label=y_label,title=title)
            mins_maxs = []
            for i, v in enumerate(y_vars):
                if v in monthly_df.columns:
                    min_max = Plot.scatter_plot(
                        fig, 'energy', v, monthly_source, colors[i], 
                        label=labels[i]
                    )
                    mins_maxs.append(min_max)
            mins_maxs = np.array(mins_maxs)
            one2one_vals = np.arange(min(mins_maxs[:,0]), max(mins_maxs[:,1]),1)
            fig.line(
                one2one_vals, one2one_vals, legend='1:1 line', color='black', 
                line_dash='dashed'
            )
            monthly_figs.append(fig)
        else:
            print('Energy balance scatter grapths missing all variables')


        #### 
        # latent energy scatter plots
        #### 
        title = 'Daily Latent Energy, Initial versus Corrected'
        unit = _get_units(['LE', 'LE_corr', 'LE_user_corr'], units)
        y_label = 'corrected ({})'.format(unit)
        x_label = 'initial ({})'.format(unit)
        fig = figure(x_axis_label=x_label, y_axis_label=y_label, title=title)
        y_vars = ['LE_corr', 'LE_user_corr']
        colors = ['red', 'darkred']
        labels = ['corr', 'user_corr']
        # add plot pairs to plot if they exist, add 1:1
        mins_maxs = []
        n_vars_fnd = 0
        for i, v in enumerate(y_vars):
            if v in df.columns and not df[v].isna().all():
                n_vars_fnd += 1
                min_max = Plot.scatter_plot(
                    fig, 'LE', v, daily_source, colors[i], label=labels[i]
                )
                mins_maxs.append(min_max)
        if n_vars_fnd > 0:
            # add scaled one to one line
            mins_maxs = np.array(mins_maxs)
            one2one_vals = np.arange(min(mins_maxs[:,0]), max(mins_maxs[:,1]),1)
            fig.line(
                one2one_vals, one2one_vals, legend='1:1 line', color='black', 
                line_dash='dashed'
            )
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Latent Energy, Initial versus Corrected'
            fig = figure(x_axis_label=x_label, y_axis_label=y_label,title=title)
            mins_maxs = []
            for i, v in enumerate(y_vars):
                if v in monthly_df.columns:
                    min_max = Plot.scatter_plot(
                        fig, 'LE', v, monthly_source, colors[i], 
                        label=labels[i]
                    )
                    mins_maxs.append(min_max)
            mins_maxs = np.array(mins_maxs)
            one2one_vals = np.arange(min(mins_maxs[:,0]), max(mins_maxs[:,1]),1)
            fig.line(
                one2one_vals, one2one_vals, legend='1:1 line', color='black', 
                line_dash='dashed'
            )
            monthly_figs.append(fig)
        else:
            print('Latent energy scatter grapths missing all variables')

        #### 
        # ET scatter plots
        #### 
        title = 'Daily Evapotranspiration, Initial versus Corrected'
        unit = _get_units(['et', 'et_corr', 'et_user_corr'], units)
        y_label = 'corrected ({})'.format(unit)
        x_label = 'initial ({})'.format(unit)
        fig = figure(x_axis_label=x_label, y_axis_label=y_label, title=title)
        y_vars = ['et_corr', 'et_user_corr']
        colors = ['red', 'darkred']
        labels = ['corr', 'user_corr']
        # add plot pairs to plot if they exist, add 1:1
        mins_maxs = []
        n_vars_fnd = 0
        for i, v in enumerate(y_vars):
            if v in df.columns and not df[v].isna().all():
                n_vars_fnd += 1
                min_max = Plot.scatter_plot(
                    fig, 'et', v, daily_source, colors[i], label=labels[i]
                )
                mins_maxs.append(min_max)
        if n_vars_fnd > 0:
            # add scaled one to one line
            mins_maxs = np.array(mins_maxs)
            one2one_vals = np.arange(min(mins_maxs[:,0]), max(mins_maxs[:,1]),1)
            fig.line(
                one2one_vals, one2one_vals, legend='1:1 line', color='black', 
                line_dash='dashed'
            )
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Evapotranspiration, Initial versus Corrected'
            fig = figure(x_axis_label=x_label, y_axis_label=y_label,title=title)
            mins_maxs = []
            for i, v in enumerate(y_vars):
                if v in monthly_df.columns:
                    min_max = Plot.scatter_plot(
                        fig, 'et', v, monthly_source, colors[i], 
                        label=labels[i]
                    )
                    mins_maxs.append(min_max)
            mins_maxs = np.array(mins_maxs)
            one2one_vals = np.arange(min(mins_maxs[:,0]), max(mins_maxs[:,1]),1)
            fig.line(
                one2one_vals, one2one_vals, legend='1:1 line', color='black', 
                line_dash='dashed'
            )
            monthly_figs.append(fig)
        else:
            print('Evapotranspiration scatter grapths missing all variables')


        #### 
        # multiple soil heat flux sensor time series plots
        #### 
        # keep user names for these in hover 
        g_re = re.compile('g_\d+|G')
        g_vars = [
            v for v in variables if g_re.match(v) and v in df.columns
        ]
        num_lines = len(g_vars)
        if num_lines > 0:
            rename_dict = {k:variables[k] for k in g_vars}
            tmp_df = df[g_vars].rename(columns=rename_dict)
            tmp_source = ColumnDataSource(tmp_df)
            plt_vars = list(rename_dict.values())
            colors = Viridis256[0:-1:int(256/num_lines)]
            title = 'Daily Soil Moisture (Multiple Sensors)'
            x_label = 'date'
            y_label = _get_units(g_vars, units)
            fig = figure(x_axis_label=x_label,y_axis_label=y_label,title=title)
            fig = Plot.add_lines(
                fig, tmp_df, plt_vars, colors, x_label, tmp_source, 
                labels=plt_vars
            )
            if fig is not None:
                daily_figs.append(fig)
                # same for monthly fig
                tmp_df = monthly_df[g_vars].rename(columns=rename_dict)
                tmp_source = ColumnDataSource(tmp_df)
                title = 'Monthly Soil Moisture (Multiple Sensors)'
                fig = figure(
                    x_axis_label=x_label, y_axis_label=y_label,title=title
                )
                fig = Plot.add_lines(
                    fig, tmp_df, plt_vars, colors, x_label, tmp_source,
                    labels=plt_vars
                )
                monthly_figs.append(fig)
            # do not print warning if missing multiple soil moisture recordings


        #### 
        # multiple soil moisture time series plots
        #### 
        # keep user names for these in hover 
        theta_re = re.compile('theta_[\d+|mean]')
        theta_vars = [
            v for v in variables if theta_re.match(v) and v in df.columns
        ]
        num_lines = len(theta_vars)
        if num_lines > 0:
            rename_dict = {k:variables[k] for k in theta_vars}
            tmp_df = df[theta_vars].rename(columns=rename_dict)
            tmp_source = ColumnDataSource(tmp_df)
            plt_vars = list(rename_dict.values())
            colors = Viridis256[0:-1:int(256/num_lines)]
            title = 'Daily Soil Moisture (Multiple Sensors)'
            x_label = 'date'
            y_label = _get_units(theta_vars, units)
            fig = figure(x_axis_label=x_label,y_axis_label=y_label,title=title)
            fig = Plot.add_lines(
                fig, tmp_df, plt_vars, colors, x_label, tmp_source, 
                labels=plt_vars
            )
            if fig is not None:
                daily_figs.append(fig)
                # same for monthly fig
                tmp_df = monthly_df[theta_vars].rename(columns=rename_dict)
                tmp_source = ColumnDataSource(tmp_df)
                title = 'Monthly Soil Moisture (Multiple Sensors)'
                fig = figure(
                    x_axis_label=x_label, y_axis_label=y_label,title=title
                )
                fig = Plot.add_lines(
                    fig, tmp_df, plt_vars, colors, x_label, tmp_source,
                    labels=plt_vars
                )
                monthly_figs.append(fig)
            # do not print warning if missing multiple soil moisture recordings




        ####
        # Aggregate plots and output depending on options
        ####
        figs = daily_figs + monthly_figs
        # remove None values in list 
        figs = list(filter(None, figs))
        grid = gridplot(
            figs, ncols=ncols, plot_width=450, plot_height=450, 
            sizing_mode='scale_width'
        )
        if output_type == 'show':
            show(grid)
        elif output_type == 'save':
            save(grid)

