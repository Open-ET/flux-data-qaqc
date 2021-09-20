# -*- coding: utf-8 -*-
import re
import numpy as np
import pandas as pd
from bokeh import models
from bokeh.io import reset_output
from bokeh.palettes import Viridis256
from bokeh.layouts import gridplot, column
from bokeh.plotting import figure, ColumnDataSource, output_file, show, save
from bokeh.models import HoverTool, Div, Range1d
from bokeh.models.formatters import DatetimeTickFormatter
from pathlib import Path


class Plot(object):
    """
    Container of plot routines of :mod:`fluxdataqaqc` including static methods
    that can be used to create and update interactive line and scatter plots
    from an arbitrary :obj:`pandas.DataFrame` instance.  
    
    Note: 
        The :obj:`.Data` and :obj:`.QaQc` objects both inherit all methods of
        :obj:`Plot` therefore allowing them to be easily used for custom
        interactive time series plots for data within input data (in
        :any:`fluxdataqaqc.Data.df`) and daily and monthly data in
        :attr:`fluxdataqaqc.QaQc.df` and :attr:`.QaQc.monthly_df`.

    """

    def __init__(self):
        pass

    @staticmethod
    def line_plot(fig, x, y, source, color, label=None,
            x_axis_type='date', **kwargs):
        """
        Add a single time series to a :obj:`bokeh.plotting.figure.Figure`
        object using data from a datetime indexed :obj:`pandas.DataFrame` with
        an interactive hover tool. 

        Interactive hover shows the values of all time series data and date
        that is added to the figure.

        Arguments:
            fig (:obj:`bokeh.plotting.figure.Figure`): a figure instance to add 
                the line to.
            x (str): name of the datetime index or column in the 
                :obj:`pandas.DataFrame` containing data to plot.
            y (str): name of the column in the :obj:`pandas.DataFrame` to plot.
            source (:obj:`bokeh.models.sources.ColumnDataSource`): column data 
                source created from the :obj:`pandas.DataFrame` with data to 
                plot.
            color (str): color of plot line, see Bokeh for color options.
            label (str or :obj:`None`): default :obj:`None`. Label for plot 
                legend (for ``y``).
            x_axis_type (:obj:`str` or :obj:`None`): default 'date'. If "date" 
                then the x-axis will be formatted as month-day-year. 

        Returns:
            :obj:`None`

        Example:

            To use the :meth:`Plot.line_plot` function we first need to create
            a :obj:`bokeh.models.sources.ColumnDataSource` from a
            :obj:`pandas.DataFrame`. Let's say we want to plot the monthly time
            series of corrected latent energy, starting from a config.ini file,

            >>> from fluxdataqaqc import Data, QaQc, Plot
            >>> d = Data('path/to/config.ini')
            >>> q = QaQc(d)
            >>> q.correct_data()
            
            Now the :obj:`.QaQc` should have the "LE_corr" (corrected latent
            energy) column, we can now make a
            :obj:`bokeh.models.sources.ColumnDataSource` from
            :attr:`fluxdataqaqc.QaQc.df` or
            :attr:`fluxdataqaqc.QaQc.monthly_df`,

            >>> from bokeh.plotting import ColumnDataSource, figure, show
            >>> source = ColumnDataSource(q.monthly_df)
            >>> # create the figure before using line_plot
            >>> fig = figure(x_axis_label='date', y_axis_label='Corrected LE')
            >>> Plot.line_plot(
            >>>     fig, 'date', 'LE_corr', source, color='red', line_width=3
            >>> )
            >>> show(fig)

            Notice, ``line_width`` is not an argument to :meth:`Plot.line_plot`
            but it is an acceptable keyword argument to
            :obj:`bokeh.plotting.figure.Figure` and therefore will work as
            expected.


        Note:
            This method is also available from the :obj:`.Data` and :obj:`.QaQc`
            objects.
        """
        hover_mode = kwargs.pop('mode','vline')
        if label is None:
            fig.line(x,y, source=source, color=color, **kwargs)
        else:
            #label=dict(value=label) # old label requirement bokeh < 2? 
            fig.line(x,y,source=source,color=color,legend_label=label,**kwargs)
            fig.legend.location = "top_left"
            fig.legend.click_policy="hide"
        # add Hover tool with additional tips if more than one line
        if len(fig.hover) == 0:
            if x_axis_type == 'date':
                Hover = HoverTool(
                    tooltips=[
                        (x,'@{}{}'.format(x,'{%F}')), 
                        (y,'@{}'.format(y))
                    ], 
                    formatters={'@{}'.format(x): 'datetime'},
                    mode = hover_mode 
                )
                fig.add_tools(Hover)
            else: 
                Hover = HoverTool(
                    tooltips=[
                        (x,'@{}'.format(x)), 
                        (y,'@{}'.format(y))
                    ], 
                    mode = hover_mode 
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
        Add paired time series data to an interactive Bokeh scatter plot 
        :obj:`bokeh.plotting.figure.Figure`.

        Handles missing data points (gaps) by masking out indices in ``x`` and
        ``y`` where one or both are null. The ``lsrl`` option adds the best
        fit least squares linear regression line with y-intercept through zero
        and reports the slope of the line in the figure legend. Interactive
        hover shows the values of all paired (x,y) data and date that is added
        to the figure. 
        
        Returns:
            (tuple): minimum and maximum ``x`` and minimum and maximum ``y`` values of paired data which can be used for adding a one to one line to the figure or other uses.

        Example:

            Let's say that we wanted to run the energy balance ratio closure
            correction including gap filling with reference ET * crop
            coefficient and then plot corrected ET versus the calculated ET
            from reference ET (named "et_fill" in ``flux-data-qaqc``) which is
            calculated on all days even those without gaps. Similar to
            :meth:`Plot.line_plot` we first need to create a
            :obj:`bokeh.models.sources.ColumnDataSource` from a
            :obj:`pandas.DataFrame`. 

            >>> from fluxdataqaqc import Data, QaQc
            >>> d = Data('path/to/config.ini')
            >>> q = QaQc(d)
            >>> q.correct_data()
            
            Now the :obj:`.QaQc` instance should have the "et_corr" (corrected
            ET) and "et_fill" (et calculated from reference ET and crop
            coefficient) columns, we can now make a
            :obj:`bokeh.models.sources.ColumnDataSource` from
            :attr:`fluxdataqaqc.QaQc.df` or
            :attr:`fluxdataqaqc.QaQc.monthly_df`,

            >>> from bokeh.plotting import ColumnDataSource, figure, show
            >>> df = q.df
            >>> source = ColumnDataSource(df)
            >>> fig = figure(
            >>>     x_axis_label='ET, corrected', y_axis_label='ET, gap fill'
            >>> )
            >>> # note, we are calling this plot method from a QaQc instace
            >>> q.scatter_plot(
            >>>     fig, 'ET_corr', 'ET_fill', source, 'red', label='lslr'
            >>> )
            >>> show(fig)

            The ``label`` keyword argument will be used in the legend and since
            the least squares linear regression line between x and y is being
            calculated the slope of the line will also be printed in the legend.
            In this case, if the slope of the regression line is 0.94 then the 
            legend will read "lslr, slope=0.94".

        Note:
            Extra keyword arguments (accepted by
            :obj:`bokeh.plotting.figure.Figure`) will be passed to the scatter
            plot but not to the least squares regression line plot.

        Note:
            This method is also available from the :obj:`.Data` and :obj:`.QaQc`
            objects.
        """
        name = '{}_vs_{}'.format(x,y)
        # remove pairs where one or both are nans, for LSRL and min/max-1:1
        
        if not y in source.data:
            print(f'WARNING: {y} not in ColumnDataSource, cannot plot.')
            return

        if not x in source.data:
            print(f'WARNING: {x} not in ColumnDataSource, cannot plot.')
            return

        xd = source.data[x].astype(float)
        yd = source.data[y].astype(float)
        mask = (~np.isnan(xd) & ~np.isnan(yd))
        if not mask.any():
            print(f'WARNING: cannot plot {name} because no paired data')
            return

        xd = xd[mask]
        yd = yd[mask]

        # allow modifications of styles, fallback
        in_a = kwargs.pop('fill_alpha', 0.2)
        size = kwargs.pop('size', 10)
        # ordinary least squares linear regression slope through zero
        if lsrl:
            #source.data.update(ColumnDataSource({x: xd, y: yd}).data)
            # plot scatter and linear regression line slope for y=mx
            m = np.linalg.lstsq(xd.reshape(-1,1),yd,rcond=None)[0][0]
            fig.scatter(
                x, y, source=source, color=color, line_width=1, fill_alpha=in_a,
                legend_label='{lab}, slope={s:.2f}'.format(lab=label, s=m), 
                name=name, size=size, **kwargs
            )
            fig.line(xd, m * xd, color=color)
        else:
            fig.scatter(
                x, y, source=source, color=color, line_width=1, fill_alpha=in_a,
                legend_label=label, name=name, size=size, **kwargs
            )
        if len(fig.hover) == 0:
            if x == date_name:
                Hover = HoverTool(
                    tooltips=[
                        (date_name,'@{}{}'.format(date_name,'{%F}')),
                        (y,'@{}'.format(y))
                    ],formatters={f'@{date_name}': 'datetime'}
                )
            else:
                Hover = HoverTool(
                    tooltips=[
                        (date_name,'@{}{}'.format(date_name,'{%F}')),
                        (x,'@{}'.format(x)), 
                        (y,'@{}'.format(y))
                    ],formatters={f'@{date_name}': 'datetime'}
                )

            fig.add_tools(Hover)
        else:
            if x == date_name:
                new_tips = [
                    (y,'@{}'.format(y))
                ]
            else:
                new_tips = [(x,'@{}'.format(x)),(y,'@{}'.format(y))]
            for t in new_tips:
                if not t in fig.hover[0].tooltips:
                    fig.hover[0].tooltips.append(t)
        fig.hover[0].names.append(name)
        fig.legend.location = "top_left"
        fig.legend.click_policy="hide"
        
        if xd.size > 0:
            return xd.min(), xd.max(), yd.min(), yd.max()

    @staticmethod
    def add_lines(fig, df, plt_vars, colors, x_name, source, 
            labels=None, **kwargs): 
        """
        Add a multiple time series to a :obj:`bokeh.plotting.figure.Figure`
        object using data from a datetime indexed :obj:`pandas.DataFrame` with
        an interactive hover tool. 

        Interactive hover shows the values of all time series data and date
        that is added to the figure.

        Arguments:
            df (:obj:`pandas.DataFrame`): :obj:`pandas.DataFrame` containing 
                time series data.
            plt_vars (list): list of data columns in ``df`` to plot.
            colors (list): list of line colors for variables in ``plt_vars``.
            x_name (str): name of the x-axis variable, e.g. the datetime index,
                in the :obj:`pandas.DataFrame` (``df``) containing data to plot.
            source (:obj:`bokeh.models.sources.ColumnDataSource`): column data 
                source created from the :obj:`pandas.DataFrame` with data to 
                plot, i.e. ``df``.
            labels (:obj:`list` or :obj:`None`): default :obj:`None`. Labels for
                each plot variable in ``plt_vars``. 

        Returns:
            ret (:obj:`None` or :obj:`bokeh.plotting.figure.Figure`): if none of the variables in ``plt_vars`` are found in ``df`` then return :obj:`None` otherwise returns the updated figure. 

        Example:

            Similar to :meth:`Plot.line_plot` we first need to create a
            :obj:`bokeh.models.sources.ColumnDataSource` from a
            :obj:`pandas.DataFrame`. This example shows how to plot two
            variables, daily corrected latent energy and sensible heat on the
            same plot. 

            >>> from fluxdataqaqc import Data, QaQc, Plot
            >>> d = Data('path/to/config.ini')
            >>> q = QaQc(d)
            >>> q.correct_data()
            
            Now the :obj:`.QaQc` instance should have the "LE_corr" (corrected
            latent energy) and "H_corr" (corrected sensible heat) columns, we
            can now make a :obj:`bokeh.models.sources.ColumnDataSource` from
            :attr:`fluxdataqaqc.QaQc.df` or
            :attr:`fluxdataqaqc.QaQc.monthly_df`,

            >>> from bokeh.plotting import ColumnDataSource, figure, show
            >>> df = q.df
            >>> plt_vars = ['LE_corr', 'H_corr']
            >>> colors = ['blue', 'red']
            >>> labels = ['LE', 'H']
            >>> source = ColumnDataSource(df)
            >>> fig = figure(
            >>>     x_axis_label='date', y_axis_label='Corrected Turbulent Flux'
            >>> )
            >>> Plot.add_lines(
            >>>     fig, df, plt_vars, colors, 'date', source, labels=labels
            >>> )
            >>> show(fig)

        Note:
            This method is also available from the :obj:`.Data` and :obj:`.QaQc`
            objects.

        """
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

    def _plot(self, FluxObj, ncols=1, output_type='save', out_file=None, 
            suptitle='', plot_width=1000, plot_height=450, 
            sizing_mode='scale_both', merge_tools=False, link_x=True, **kwargs): 
        """ 
        Private routine for aggregated validation plots that are used by
        the :meth:`.QaQc.plot` and :meth:`.Data.plot` methods.
        """
        # get daily and monthly time series with internal names, get units
        monthly = False
        if hasattr(FluxObj, 'monthly_df'):
            # will run correction as of now if it is a QaQc
            monthly = True
            monthly_df = FluxObj.monthly_df.rename(columns=FluxObj.inv_map) 
            # avoid plotting single point- errors out bokeh datetime axis, etc.
            for c in monthly_df.columns:
                if monthly_df[c].notna().sum() <= 1:
                    monthly_df.drop(c, axis=1, inplace=True)
            monthly_source = ColumnDataSource(monthly_df)

        # so that the correction is run, may change this
        FluxObj.df.head(); # if Data, need to access to calc vp/vpd 
        df = FluxObj.df.rename(columns=FluxObj.inv_map) 
        variables = FluxObj.variables
        units = FluxObj.units 
        # bokeh column sources for tooltips
        daily_source=ColumnDataSource(df)
        # for aggregating plots
        daily_line = []
        daily_scatter = []
        monthly_line = []
        monthly_scatter = []

        if output_type == 'save':
            output_file(out_file)

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


        # run through each plot, daily then monthly versions
        #### 
        # energy balance time series plots
        #### 
        plt_vars = ['LE', 'H', 'Rn', 'G']
        colors = ['blue', 'red', 'black', 'green']
        title = 'Daily Surface Energy Balance Components'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height, name='energy_balance_daily'
        )
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=plt_vars
        )
        if fig is not None:
            daily_line.append(fig)
        else:
            print(
                'Energy balance components time series grapths missing all '
                'variables'
            )
        if fig is not None and monthly:
            # same for monthly fig
            title = 'Monthly Surface Energy Balance Components'
            fig = figure(x_axis_label=x_label, y_axis_label=y_label,title=title,
                width=plot_width, height=plot_height, 
                name='energy_balance_monthly'
            )
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=plt_vars
            )
            monthly_line.append(fig)

        #### 
        # incoming shortwave and ASCE potential clear sky time series plots
        #### 
        plt_vars = ['sw_in', 'rso']
        # only plot if we have both
        if set(plt_vars).issubset(df.columns):
            labels = ['Station Rs', 'ASCE Rso']
            colors = ['black', 'red']
            title =\
                'Daily Incoming Shortwave (Rs) and ASCE Clear Sky Shortwave '+\
                'Radiation (Rso)'
            x_label = 'date'
            y_label = _get_units(plt_vars, units)
            fig = figure(x_axis_label=x_label, y_axis_label=y_label, 
                title=title, width=plot_width, height=plot_height, 
                name='Rs_daily'
            )
            fig = Plot.add_lines(
                fig, df, plt_vars, colors, x_label, daily_source, labels=labels
            )
            if fig is not None:
                daily_line.append(fig)
                ## same for monthly fig (removed for now)
                #title='Monthly Incoming Shortwave and ASCE Potential Radiation'
                #fig = figure(
                #    x_axis_label=x_label,y_axis_label=y_label,title=title,
                #    width=plot_width, height=plot_height
                #)
                #fig = Plot.add_lines(
                #    fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                #    labels=labels
                #)
                #monthly_line.append(fig)
        else:
            print(
                'Shortwave and potential clear sky radiation time series '
                'grapths missing all variables'
            )

        #### 
        # multiple soil heat flux sensor time series plots
        #### 
        # keep user names for these in hover 
        g_re = re.compile('^[gG]_[\d+mean|corr]|G$')
        g_vars = [
            v for v in variables if g_re.match(v) and v in df.columns
        ]
        num_lines = len(g_vars)
        if num_lines > 1:
            rename_dict = {k:variables[k] for k in g_vars}
            tmp_df = df[g_vars].rename(columns=rename_dict)
            tmp_source = ColumnDataSource(tmp_df)
            plt_vars = list(rename_dict.values())
            colors = Viridis256[0:-1:int(256/num_lines)]
            title = 'Daily Soil Heat Flux (Multiple Sensors)'
            x_label = 'date'
            y_label = _get_units(g_vars, units)
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                plot_width=plot_width, plot_height=plot_height, name='G_daily'
            )
            fig = Plot.add_lines(
                fig, tmp_df, plt_vars, colors, x_label, tmp_source, 
                labels=plt_vars
            )
            if fig is not None:
                daily_line.append(fig)
            if fig is not None and monthly:
                # same for monthly fig
                g_vars = [
                    v for v in variables if g_re.match(v) and v in \
                        monthly_df.columns
                ]
                num_lines = len(g_vars)
                if num_lines > 1:
                    tmp_df = monthly_df[g_vars].rename(columns=rename_dict)
                    tmp_source = ColumnDataSource(tmp_df)
                    title = 'Monthly Soil Heat Flux (Multiple Sensors)'
                    fig = figure(
                        x_axis_label=x_label, y_axis_label=y_label,title=title,
                        plot_width=plot_width, plot_height=plot_height, 
                        name='G_monthly'
                    )
                    fig = Plot.add_lines(
                        fig, tmp_df, plt_vars, colors, x_label, tmp_source,
                        labels=plt_vars
                    )
                    monthly_line.append(fig)
            # do not print warning if missing multiple soil moisture recordings

        #### 
        # radiation time series plots
        #### 
        plt_vars = ['sw_in', 'lw_in', 'sw_out', 'lw_out']
        colors = ['red', 'darkred', 'blue', 'navy']
        title = 'Daily Radiation Components'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height, name='radiation_daily'
        )
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=plt_vars
        )
        if fig is not None:
            daily_line.append(fig)
        else:
            print(
                'Radiation components time series grapths missing all variables'
            )
        if fig is not None and monthly:
            # same for monthly fig
            title = 'Monthly Radiation Components'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height, name='radiation_monthly'
            )
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=plt_vars
            )
            monthly_line.append(fig)


        #### 
        # average temperature time series plot
        #### 
        plt_vars = ['t_avg']
        colors = ['black']
        title = 'Daily Average Air Temperature'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height, name='temp_daily'
        )
        fig = Plot.add_lines(fig, df, plt_vars, colors, x_label, daily_source)
        if fig is not None:
            daily_line.append(fig)
        else:
            print(
                'Average air temperature time series grapths missing all '
                'variables'
            )
        if fig is not None and monthly:
            # same for monthly fig
            title = 'Monthly Average Air Temperature'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label,title=title,
                width=plot_width, height=plot_height, name='temp_monthly'
            )
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source
            )
            monthly_line.append(fig)

        #### 
        # vapor pressure time series plots
        #### 
        plt_vars = ['vp', 'vpd']
        colors = ['black', 'darkred']
        title = 'Daily Average Vapor Pressure and Deficit'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height, name='vap_press_daily'
        )
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=plt_vars
        )
        if fig is not None:
            daily_line.append(fig)
        else:
            print('Vapor pressure time series grapths missing all variables')
        if fig is not None and monthly:
            # same for monthly fig
            title = 'Monthly Average Vapor Pressure'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height, name='vap_press_monthly'
            )
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=plt_vars
            )
            monthly_line.append(fig)

        #### 
        # windpseed time series plot
        #### 
        plt_vars = ['ws']
        colors = ['black']
        title = 'Daily Average Windspeed'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height, name='wind_daily'
        )
        fig = Plot.add_lines(fig, df, plt_vars, colors, x_label, daily_source)
        if fig is not None:
            daily_line.append(fig)
        else:
            print('Windspeed time series grapths missing all variables')
        if fig is not None and monthly:
            # same for monthly fig
            title = 'Monthly Average Windspeed'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height, name='wind_monthly'
            )
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source
            )
            monthly_line.append(fig)

        #### 
        # precipitation time series plots
        #### 
        plt_vars = ['ppt', 'gridMET_prcp']
        labels = ['station', 'gridMET']
        colors = ['black', 'red']
        title = 'Daily Station and gridMET Precipitation'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height, name='precip_daily'
        )
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=labels
        )
        if fig is not None:
            daily_line.append(fig)
        else:
            print('Precipitation time series grapths missing all variables')
        if fig is not None and monthly:
            # same for monthly fig
            title = 'Monthly Station and gridMET Precipitation'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height, name='precip_monthly'
            )
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=labels
            )
            monthly_line.append(fig)

        #### 
        # latent energy time series plots
        #### 
        plt_vars = ['LE', 'LE_corr', 'LE_user_corr']
        colors = ['black', 'red', 'darkorange']
        title = 'Daily Average Latent Energy Flux'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height, name='LE_daily'
        )
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=plt_vars
        )
        if fig is not None:
            daily_line.append(fig)
        else:
            print('Latent energy time series grapths missing all variables')
        if fig is not None and monthly:
            # same for monthly fig
            title = 'Monthly Average Latent Energy Flux'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height, name='LE_monthly'
            )
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=plt_vars
            )
            monthly_line.append(fig)

        #### 
        # ET time series plots
        #### 
        refET = 'ETr' if 'ETrF' in df.columns else 'ETo'
        plt_vars = ['ET', 'ET_corr', 'ET_user_corr', f'gridMET_{refET}']
        labels = plt_vars[0:3] + [refET]
        colors = ['black', 'red', 'darkorange', 'blue']
        title = 'Daily Evapotranspiration'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height, name='ET_daily'
        )
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=labels
        )
        if 'ET_fill_val' in df.columns and fig is not None:
            # make gap fill values more visible
            Plot.line_plot(
                fig, 'date', 'ET_fill_val', daily_source, 'green', 
                label='ET_fill_val', line_width=3
            )

        if fig is not None:
            daily_line.append(fig)
        else:
            print(
                'Evapotranspiration time series grapths missing all variables'
            )
        if fig is not None and monthly:
            # same for monthly fig
            title = 'Monthly Evapotranspiration'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height, name='ET_monthly'
            )
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=labels
            )
            monthly_line.append(fig)

        #### 
        # number gap filled days monthly time series plot
        #### 
        if monthly and 'ET_gap' in monthly_df.columns:
            txt = ''
            if 'ET_corr' in df.columns:
                txt = ' Corrected'
            title = 'Number of Gap Filled Days in{} Monthly ET'.format(txt)
            x_label = 'date'
            y_label = 'number of gap-filled days'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height, name='ET_gaps'
            )
            x = 'date'
            y = 'ET_gap'
            color = 'black'
            Plot.line_plot(fig, x, y, monthly_source, color)
            monthly_line.append(fig)
        elif monthly:
            print('Monthly count of gap filled ET days plot missing variable')

        #### 
        # ETrF time series plots
        ####
        plt_vars = [f'{refET}F', f'{refET}F_filtered']
        colors = ['black', 'red']
        title = f'Daily Fraction of Reference ET ({refET}F)'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height, name=f'{refET}F_daily'
        )
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=plt_vars
        )
        if fig is not None:
            daily_line.append(fig)
        else:
            print(
                'Fraction of reference ET time series grapths missing all '
                'variables'
            )
        if fig is not None and monthly:
            # same for monthly fig
            title = f'Monthly Fraction of Reference ET ({refET}F)'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height, name=f'{refET}F_monthly'
            )
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=plt_vars
            )
            monthly_line.append(fig)

        #### 
        # energy balance ratio time series plots
        #### 
        plt_vars = ['ebr', 'ebr_corr', 'ebr_user_corr']
        colors = ['black', 'red', 'darkorange']
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
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height, name='EBR_daily'
        )
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=labels
        )
        if fig is not None:
            daily_line.append(fig)
        else:
            print(
                'Energy balance ratio time series grapths missing all '
                'variables'
            )
        if fig is not None and monthly:
            # same for monthly fig
            title = 'Monthly Energy Balance Ratio with Long-term Mean'
            # add mean for monthly EBRs to legend
            labels = []
            for i, v in enumerate(plt_vars):
                if v in monthly_df.columns:
                    added_text = ': {}'.format(
                        str(round(monthly_df[v].mean(),2))
                    )
                    labels.append(plt_vars[i] + added_text)
                else:
                    labels.append(None)
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height, name='EBR_monthly'
            )
            fig = Plot.add_lines(
                fig, monthly_df, plt_vars, colors, x_label, monthly_source,
                labels=labels
            )
            monthly_line.append(fig)

        #### 
        # energy balance closure scatter plots
        #### 
        title = 'Daily Energy Balance Closure, Energy Versus Flux with Slope '\
            'Through Origin'
        unit = _get_units(['LE', 'H', 'Rn', 'G'], units)
        y_label = 'LE + H ({})'.format(unit)
        x_label = 'Rn - G ({})'.format(unit)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_width, name='energy_vs_flux_daily'
        )
        y_vars = ['flux', 'flux_corr', 'flux_user_corr']
        colors = ['black', 'red', 'darkorange']
        labels = ['init', 'corr', 'user_corr']
        # add plot pairs to plot if they exist, add 1:1
        mins_maxs = []
        n_vars_fnd = 0
        for i, v in enumerate(y_vars):
            if v in df.columns and not df[v].isna().all():
                n_vars_fnd += 1
                if v == 'flux_corr' and 'energy_corr' in df.columns:
                    x_var = 'energy_corr'
                else:
                    x_var = 'energy'
                min_max = Plot.scatter_plot(
                    fig, x_var, v, daily_source, colors[i], label=labels[i]
                )
                if min_max is not None:
                    mins_maxs.append(min_max)
        if n_vars_fnd > 0:
            # add scaled one to one line
            mins_maxs = np.array(mins_maxs)
            if not pd.isna(mins_maxs).all():
                x_min = min(mins_maxs[:,0])
                x_max = max(mins_maxs[:,1])
                y_min = min(mins_maxs[:,2])
                y_max = max(mins_maxs[:,3])
                ax_min, ax_max = min([x_min,y_min]), max([x_max,y_max])
                ax_min -= 0.02*abs(ax_max-ax_min)
                ax_max += 0.02*abs(ax_max-ax_min)
                fig.x_range=Range1d(ax_min, ax_max)
                fig.y_range=Range1d(ax_min, ax_max)
                one2one_vals = np.arange(ax_min, ax_max,1)
                fig.line(
                    one2one_vals, one2one_vals, legend_label='1:1 line', 
                    color='black', line_dash='dashed'
                )
                daily_scatter.append(fig)
            if monthly:
                # same for monthly fig
                title = 'Monthly Energy Balance Closure, Energy Versus Flux '\
                    'with Slope Through Origin'
                fig = figure(
                    x_axis_label=x_label, y_axis_label=y_label, title=title,
                    width=plot_width, height=plot_width, 
                    name='energy_vs_flux_monthly'
                )
                mins_maxs = []
                for i, v in enumerate(y_vars):
                    if v in monthly_df.columns:
                        min_max = Plot.scatter_plot(
                            fig, 'energy', v, monthly_source, colors[i], 
                            label=labels[i]
                        )
                        if min_max is not None:
                            mins_maxs.append(min_max)
                mins_maxs = np.array(mins_maxs)
                # check if not all pairs are empty, if not plot 1:1
                if not pd.isna(mins_maxs).all():
                    x_min = min(mins_maxs[:,0])
                    x_max = max(mins_maxs[:,1])
                    y_min = min(mins_maxs[:,2])
                    y_max = max(mins_maxs[:,3])
                    ax_min, ax_max = min([x_min,y_min]), max([x_max,y_max])
                    ax_min -= 0.02*abs(ax_max-ax_min)
                    ax_max += 0.02*abs(ax_max-ax_min)
                    fig.x_range=Range1d(ax_min, ax_max)
                    fig.y_range=Range1d(ax_min, ax_max)
                    one2one_vals = np.arange(ax_min, ax_max,1)
                    fig.line(
                        one2one_vals, one2one_vals, legend_label='1:1 line', 
                        color='black', line_dash='dashed'
                    )
                    monthly_scatter.append(fig)
        else:
            print('Energy balance scatter grapths missing all variables')


        #### 
        # latent energy scatter plots
        #### 
        title = 'Daily Latent Energy, Initial Versus Corrected'
        unit = _get_units(['LE', 'LE_corr', 'LE_user_corr'], units)
        y_label = 'corrected ({})'.format(unit)
        x_label = 'initial ({})'.format(unit)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_width, name='LE_scatter_daily'
        )
        y_vars = ['LE_corr', 'LE_user_corr']
        colors = ['red', 'darkorange']
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
            if not pd.isna(mins_maxs).all():
                x_min = min(mins_maxs[:,0])
                x_max = max(mins_maxs[:,1])
                y_min = min(mins_maxs[:,2])
                y_max = max(mins_maxs[:,3])
                ax_min, ax_max = min([x_min,y_min]), max([x_max,y_max])
                ax_min -= 0.02*abs(ax_max-ax_min)
                ax_max += 0.02*abs(ax_max-ax_min)
                fig.x_range=Range1d(ax_min, ax_max)
                fig.y_range=Range1d(ax_min, ax_max)
                one2one_vals = np.arange(ax_min, ax_max,1)
                fig.line(
                    one2one_vals, one2one_vals, legend_label='1:1 line', 
                    color='black', line_dash='dashed'
                )
                daily_scatter.append(fig)
            if monthly:
                # same for monthly fig
                title = 'Monthly Latent Energy, Initial Versus Corrected'
                fig = figure(
                    x_axis_label=x_label, y_axis_label=y_label, title=title,
                    width=plot_width, height=plot_width, 
                    name='LE_scatter_monthly'
                )
                mins_maxs = []
                for i, v in enumerate(y_vars):
                    if v in monthly_df.columns:
                        min_max = Plot.scatter_plot(
                            fig, 'LE', v, monthly_source, colors[i], 
                            label=labels[i]
                        )
                        if min_max is not None:
                            mins_maxs.append(min_max)
                mins_maxs = np.array(mins_maxs)
                # check if not all pairs are empty, if not plot 1:1
                if not pd.isna(mins_maxs).all():
                    x_min = min(mins_maxs[:,0])
                    x_max = max(mins_maxs[:,1])
                    y_min = min(mins_maxs[:,2])
                    y_max = max(mins_maxs[:,3])
                    ax_min, ax_max = min([x_min,y_min]), max([x_max,y_max])
                    ax_min -= 0.02*abs(ax_max-ax_min)
                    ax_max += 0.02*abs(ax_max-ax_min)
                    fig.x_range=Range1d(ax_min, ax_max)
                    fig.y_range=Range1d(ax_min, ax_max)
                    one2one_vals = np.arange(ax_min, ax_max,1)
                    fig.line(
                        one2one_vals, one2one_vals, legend_label='1:1 line', 
                        color='black', line_dash='dashed'
                    )
                    monthly_scatter.append(fig)
        else:
            print('Latent energy scatter grapths missing all variables')

        #### 
        # ET scatter plots
        #### 
        title = 'Daily Evapotranspiration, Initial Versus Corrected'
        unit = _get_units(['ET', 'ET_corr', 'ET_user_corr'], units)
        y_label = 'corrected ({})'.format(unit)
        x_label = 'initial ({})'.format(unit)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_width, name='ET_scatter_daily'
        )
        y_vars = ['ET_corr', 'ET_user_corr']
        colors = ['red', 'darkorange']
        labels = ['corr', 'user_corr']
        # add plot pairs to plot if they exist, add 1:1
        mins_maxs = []
        n_vars_fnd = 0
        for i, v in enumerate(y_vars):
            if v in df.columns and not df[v].isna().all():
                n_vars_fnd += 1
                min_max = Plot.scatter_plot(
                    fig, 'ET', v, daily_source, colors[i], label=labels[i]
                )
                mins_maxs.append(min_max)
        if n_vars_fnd > 0:
            # add scaled one to one line
            mins_maxs = np.array(mins_maxs)
            x_min = min(mins_maxs[:,0])
            x_max = max(mins_maxs[:,1])
            y_min = min(mins_maxs[:,2])
            y_max = max(mins_maxs[:,3])
            ax_min, ax_max = min([x_min,y_min]), max([x_max,y_max])
            ax_min -= 0.02*abs(ax_max-ax_min)
            ax_max += 0.02*abs(ax_max-ax_min)
            fig.x_range=Range1d(ax_min, ax_max)
            fig.y_range=Range1d(ax_min, ax_max)
            one2one_vals = np.arange(ax_min, ax_max,1)
            fig.line(
                one2one_vals, one2one_vals, legend_label='1:1 line', 
                color='black', line_dash='dashed'
            )
            daily_scatter.append(fig)
            if monthly:
                # same for monthly fig
                title = 'Monthly Evapotranspiration, Initial Versus Corrected'
                fig = figure(
                    x_axis_label=x_label, y_axis_label=y_label, title=title,
                    width=plot_width, height=plot_width, 
                    name='ET_scatter_monthly'
                )
                mins_maxs = []
                for i, v in enumerate(y_vars):
                    if v in monthly_df.columns:
                        min_max = Plot.scatter_plot(
                            fig, 'ET', v, monthly_source, colors[i], 
                            label=labels[i]
                        )
                        mins_maxs.append(min_max)
                mins_maxs = np.array(mins_maxs)
                # check if not all pairs are empty, if not plot 1:1
                if not pd.isna(mins_maxs).all():
                    x_min = min(mins_maxs[:,0])
                    x_max = max(mins_maxs[:,1])
                    y_min = min(mins_maxs[:,2])
                    y_max = max(mins_maxs[:,3])
                    ax_min, ax_max = min([x_min,y_min]), max([x_max,y_max])
                    ax_min -= 0.02*abs(ax_max-ax_min)
                    ax_max += 0.02*abs(ax_max-ax_min)
                    fig.x_range=Range1d(ax_min, ax_max)
                    fig.y_range=Range1d(ax_min, ax_max)
                    one2one_vals = np.arange(ax_min, ax_max,1)
                    fig.line(
                        one2one_vals, one2one_vals, legend_label='1:1 line', 
                        color='black', line_dash='dashed'
                    )
                    monthly_scatter.append(fig)
        else:
            print('Evapotranspiration scatter grapths missing all variables')

        #### 
        # multiple soil moisture time series plots
        #### 
        # keep user names for these in hover 
        theta_re = re.compile('theta_[\d+|mean]')
        theta_vars = [
            v for v in variables if theta_re.match(v) and v in df.columns
        ]
        num_lines = len(theta_vars)
        if num_lines > 0 and not df[theta_vars].isna().all().all():
            rename_dict = {k:variables[k] for k in theta_vars}
            tmp_df = df[theta_vars].rename(columns=rename_dict)
            tmp_source = ColumnDataSource(tmp_df)
            plt_vars = list(rename_dict.values())
            colors = Viridis256[0:-1:int(256/num_lines)]
            title = 'Daily Soil Moisture (Multiple Sensors)'
            x_label = 'date'
            y_label = _get_units(theta_vars, units)
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                plot_width=plot_width, plot_height=plot_height, 
                name='theta_daily'
            )
            fig = Plot.add_lines(
                fig, tmp_df, plt_vars, colors, x_label, tmp_source, 
                labels=plt_vars
            )
            if fig is not None:
                daily_line.append(fig)
            theta_vars = [
                v for v in variables if theta_re.match(v) and v in\
                    df.columns
            ]
            if fig is not None and monthly and len(theta_vars) > 0:
                # same for monthly fig
                tmp_df = monthly_df[theta_vars].rename(columns=rename_dict)
                tmp_source = ColumnDataSource(tmp_df)
                title = 'Monthly Soil Moisture (Multiple Sensors)'
                fig = figure(
                    x_axis_label=x_label, y_axis_label=y_label, title=title,
                    plot_width=plot_width, plot_height=plot_height,
                    name='theta_monthly'
                )
                fig = Plot.add_lines(
                    fig, tmp_df, plt_vars, colors, x_label, tmp_source,
                    labels=plt_vars
                )
                monthly_line.append(fig)
            # do not print warning if missing multiple soil moisture recordings


        # Aggregate plots and output depending on options
        # remove None values in different figure groups 
        daily_line = list(filter(None, daily_line))
        daily_scatter = list(filter(None, daily_scatter))
        monthly_line = list(filter(None, monthly_line))
        monthly_scatter = list(filter(None, monthly_scatter))
        # link axes for time series plots
        if link_x:
            for each in daily_line:
                each.x_range = daily_line[0].x_range
            for each in monthly_line:
                each.x_range = monthly_line[0].x_range
        figs = daily_line + daily_scatter + monthly_line + monthly_scatter
        grid = gridplot(
            figs, ncols=ncols, plot_width=None, plot_height=None, 
            sizing_mode=sizing_mode, merge_tools=merge_tools, **kwargs
        )
        if output_type == 'show':
            show(column(Div(text=suptitle),grid))
        elif output_type == 'notebook':
            from bokeh.io import output_notebook
            output_notebook()
            show(column(Div(text=suptitle),grid))
        elif output_type == 'save':
            save(column(Div(text=suptitle),grid))
        elif output_type == 'return_figs':
            return figs
        elif output_type == 'return_grid':
            return grid

        reset_output()

