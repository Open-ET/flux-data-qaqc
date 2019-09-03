# -*- coding: utf-8 -*-
import re
import numpy as np
from bokeh import models
from bokeh.palettes import Viridis256
from bokeh.layouts import gridplot, column
from bokeh.plotting import figure, ColumnDataSource, output_file, show, save
from bokeh.models import HoverTool, Div
from bokeh.models.formatters import DatetimeTickFormatter
from pathlib import Path

class Plot(object):
    """
    Container of plot routines of :mod:`fluxdataqaqc` including static methods
    that can be used to create and update interactive line and scatter plots
    from a :obj:`pandas.DataFrame` instance.  
    
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
        Add paired time series data to an interactive Bokeh scatter plot 
        :obj:`bokeh.plotting.figure.Figure`.

        Handles missing data points (gaps) by masking out indices in ``x`` and
        ``y`` where one or both are null. The ``lsrl`` option adds the best
        fit least squares linear regression line with y-intercept through zero
        and reports the slope of the line in the figure legend. Interactive
        hover shows the values of all paired (x,y) data and date that is added
        to the figure. 
        
        Returns:
            (tuple): minimum and maximum ``x`` value of paired data which can be used for adding a one to one line to the figure later.

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
            >>>     fig, 'et_corr', 'et_fill', source, 'red', label='lslr'
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
        xd = source.data[x].astype(float)
        yd = source.data[y].astype(float)
        mask = (~np.isnan(xd) & ~np.isnan(yd))
        if not mask.any():
            print('WARNING: cannot plot {} because no paired data'.format(name))
            return
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
        
        if xd.size > 0:
            return xd.min(), xd.max()

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

    def _plot(self, QaQc, ncols=1, output_type='save', out_file=None, 
            suptitle=None, plot_width=1000, plot_height=450, 
            sizing_mode='scale_both', merge_tools=False, **kwargs): 
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
        if out_file is None and output_type == 'save':
            out_file = Path(QaQc.out_dir)/'{}_plots.html'.format(QaQc.site_id)
            out_dir = out_file.parent
            if not out_dir.is_dir():
                out_dir.mkdir(parents=True, exist_ok=True)
        # to allow making any subdir that does not yet exist
        # if out_file is to a non-existent directory
        elif out_file is not None and output_type == 'save':
            out_dir = Path(out_file).parent
            if not out_dir.is_dir():
                out_dir.mkdir(parents=True, exist_ok=True)

        if output_type == 'save':
            output_file(out_file)

        # bokeh column sources for tooltips
        daily_source=ColumnDataSource(df)
        monthly_source = ColumnDataSource(monthly_df)
        # for aggregating plots
        daily_figs = []
        monthly_figs = []

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
            width=plot_width, height=plot_height
        )
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=plt_vars
        )
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Surface Energy Balance Components'
            fig = figure(x_axis_label=x_label, y_axis_label=y_label,title=title,
                width=plot_width, height=plot_height
            )
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
        labels = ['Station Rs', 'ASCE Rso']
        colors = ['black', 'red']
        title =\
            'Daily Incoming Shortwave (Rs) and ASCE Clear Sky Shortwave '+\
            'Radiation (Rso)'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height
        )
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=labels
        )
        if fig is not None:
            daily_figs.append(fig)
            ## same for monthly fig (removed for now)
            #title = 'Monthly Incoming Shortwave and ASCE Potential Radiation'
            #fig = figure(x_axis_label=x_label,y_axis_label=y_label,title=title,
            #    width=plot_width, height=plot_height
            #)
            #fig = Plot.add_lines(
            #    fig, monthly_df, plt_vars, colors, x_label, monthly_source,
            #    labels=labels
            #)
            #monthly_figs.append(fig)
        else:
            print(
                'Shortwave and potential clear sky shortwave time series '
                'grapths missing all variables'
            )

        #### 
        # multiple soil heat flux sensor time series plots
        #### 
        # keep user names for these in hover 
        g_re = re.compile('g_\d+|G')
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
                plot_width=plot_width, plot_height=plot_height
            )
            fig = Plot.add_lines(
                fig, tmp_df, plt_vars, colors, x_label, tmp_source, 
                labels=plt_vars
            )
            if fig is not None:
                daily_figs.append(fig)
                # same for monthly fig
                tmp_df = monthly_df[g_vars].rename(columns=rename_dict)
                tmp_source = ColumnDataSource(tmp_df)
                title = 'Monthly Soil Heat Flux (Multiple Sensors)'
                fig = figure(
                    x_axis_label=x_label, y_axis_label=y_label,title=title,
                    plot_width=plot_width, plot_height=plot_height
                )
                fig = Plot.add_lines(
                    fig, tmp_df, plt_vars, colors, x_label, tmp_source,
                    labels=plt_vars
                )
                monthly_figs.append(fig)
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
            width=plot_width, height=plot_height
        )
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=plt_vars
        )
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Radiation Components'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height
            )
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
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height
        )
        fig = Plot.add_lines(fig, df, plt_vars, colors, x_label, daily_source)
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Average Air Temperature'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label,title=title,
                width=plot_width, height=plot_height
            )
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
        #### TODO: make vpd on secondary y-axis 
        plt_vars = ['vp', 'vpd']
        colors = ['blue', 'black']
        title = 'Daily Average Vapor Pressure and Deficit'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height
        )
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=plt_vars
        )
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Average Vapor Pressure'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height
            )
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
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height
        )
        fig = Plot.add_lines(fig, df, plt_vars, colors, x_label, daily_source)
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Average Windspeed'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height
            )
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
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height
        )
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=labels
        )
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Station and gridMET precipitation'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height
            )
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
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height
        )
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=plt_vars
        )
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Average Latent Energy Flux'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height
            )
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
        plt_vars = ['et', 'et_corr', 'et_user_corr', 'gridMET_etr_mm']
        labels = plt_vars[0:3] + ['etr']
        colors = ['black', 'red', 'darkred', 'blue']
        title = 'Daily Evapotranspiration'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height
        )
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=labels
        )
        if 'et_fill_val' in df.columns:
            # make gap fill values more visible
            Plot.line_plot(
                fig, 'date', 'et_fill_val', daily_source, 'green', 
                label='et_fill_val', line_width=3
            )

        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Evapotranspiration'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height
            )
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
        # number gap filled days monthly time series plot
        #### 

        if 'et_gap' in monthly_df.columns:
            title = 'Number of Gap Filled Days in Corrected Monthly ET'
            x_label = 'date'
            y_label = 'number of gap-filled days'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height
            )
            x = 'date'
            y = 'et_gap'
            color = 'black'
            Plot.line_plot(fig, x, y, monthly_source, color)
            monthly_figs.append(fig)
        else:
            print('Monthly count of gap filled ET days plot missing variable')

        #### 
        # Kc time series plots
        #### 
        plt_vars = ['Kc', 'Kc_7day_mean']
        colors = ['black', 'red']
        title = 'Daily Crop Coefficient (Kc)'
        x_label = 'date'
        y_label = _get_units(plt_vars, units)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height
        )
        fig = Plot.add_lines(
            fig, df, plt_vars, colors, x_label, daily_source, labels=plt_vars
        )
        if fig is not None:
            daily_figs.append(fig)
            # same for monthly fig
            title = 'Monthly Crop Coefficient (Kc)'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height
            )
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
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_height
        )
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
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_height
            )
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
        title = 'Daily Energy Balance Closure, Energy Versus Flux with Slope '\
            'Through Origin'
        unit = _get_units(['LE', 'H', 'Rn', 'G'], units)
        y_label = 'LE + H ({})'.format(unit)
        x_label = 'Rn - G ({})'.format(unit)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_width
        )
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
                if min_max is not None:
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
            title = 'Monthly Energy Balance Closure, Energy Versus Flux '\
                'with Slope Through Origin'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_width
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
        title = 'Daily Latent Energy, Initial Versus Corrected'
        unit = _get_units(['LE', 'LE_corr', 'LE_user_corr'], units)
        y_label = 'corrected ({})'.format(unit)
        x_label = 'initial ({})'.format(unit)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_width
        )
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
            title = 'Monthly Latent Energy, Initial Versus Corrected'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_width
            )
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
        title = 'Daily Evapotranspiration, Initial Versus Corrected'
        unit = _get_units(['et', 'et_corr', 'et_user_corr'], units)
        y_label = 'corrected ({})'.format(unit)
        x_label = 'initial ({})'.format(unit)
        fig = figure(
            x_axis_label=x_label, y_axis_label=y_label, title=title,
            width=plot_width, height=plot_width
        )
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
            title = 'Monthly Evapotranspiration, Initial Versus Corrected'
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                width=plot_width, height=plot_width
            )
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
            fig = figure(
                x_axis_label=x_label, y_axis_label=y_label, title=title,
                plot_width=plot_width, plot_height=plot_height
            )
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
                    x_axis_label=x_label, y_axis_label=y_label, title=title,
                    plot_width=plot_width, plot_height=plot_height
                )
                fig = Plot.add_lines(
                    fig, tmp_df, plt_vars, colors, x_label, tmp_source,
                    labels=plt_vars
                )
                monthly_figs.append(fig)
            # do not print warning if missing multiple soil moisture recordings


        # Aggregate plots and output depending on options
        figs = daily_figs + monthly_figs
        # remove None values in list 
        figs = list(filter(None, figs))
        grid = gridplot(
            figs, ncols=ncols, plot_width=None, plot_height=None, 
            sizing_mode=sizing_mode, merge_tools=merge_tools, **kwargs
        )
        if output_type == 'show':
            show(column(Div(text=suptitle),grid))
        elif output_type == 'save':
            save(column(Div(text=suptitle),grid))

