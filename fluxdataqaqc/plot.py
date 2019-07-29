# -*- coding: utf-8 -*-
"""
Plot fluxnet data using bokeh packages
"""

from .qaqc import QaQc
from bokeh.palettes import Viridis256
from bokeh.layouts import gridplot
from bokeh.models import Label, ColumnDataSource, HoverTool
from bokeh.plotting import figure, output_file, save
import numpy as np
from math import ceil


class Plot(object):
    """
    Takes a QaQc object from qaqc.py, determines what variables are available, and makes bokeh plots for the different
    climate variables.

    """

    def __init__(self, qaqc=None):

        if isinstance(qaqc, QaQc):
            self._monthly_df = qaqc.monthly_df
            self._df = qaqc.df
            self.site_id = qaqc.site_id
            self.out_dir = qaqc.out_dir
            self.variables = qaqc.variables
            self.inv_map = qaqc.inv_map
            self.provided_vars = self._inventory_variables()

        elif qaqc is not None:
            raise TypeError("Must assign a fluxdataqaqc.qaqc.QaQc object")
        else:
            self._df = None

        self.plots_created = False
        self.plot_file = None

    def _inventory_variables(self):
        """
        Checks to see what variables are present from QaQc object by iterating through provided dictionary

        Returns:
            provided_vars (list): list of variable names that are present

        """

        provided_vars = []

        for k,v in self.variables.items():
            if v != 'na':
                provided_vars.append(k) 
            if k == 'G' and v == 'G_mean':
                provided_vars.append('G_mean')

        return provided_vars

    def _generate_plot_features(self, code, usage=''):
        """
            Generates plot features depending on what code is passed, currently a 1 at the end (ex. code 91) indicates
            that a variable is missing,

            Parameters:
                code (int): code passed by main script that indicates what type of data has been passed
                usage (str): additional text indicating why this plot is being created, it may be blank.

            Returns:
                units (str): units for plot of variable(s)
                title (str): title of plot
                var_one_name (str): first variable name
                var_one_color (str): string of color code to use for plotting variable one
                var_two_name (str): second variable name
                var_two_color (str): string of color code to use for plotting variable two
                var_three_name (str): third variable name
                var_three_color :(str) string of color code to use for plotting variable three
                var_four_name (str): ourth variable name
                var_four_color (str): string of color code to use for plotting variable four
        """

        if code == 1:  # First graph, four variables
            var_one_name = 'Rn'
            var_one_color = 'black'
            var_two_name = 'LE'
            var_two_color = 'blue'
            var_three_name = 'H'
            var_three_color = 'red'
            var_four_name = 'G'
            var_four_color = 'green'
            units = 'w/m2'
            title = usage + ' Surface Balance Eq. Components'

        elif code == 2:  # Second graph, four variables
            var_one_name = 'SW_in'
            var_one_color = 'red'
            var_two_name = 'LW_in'
            var_two_color = 'darkred'
            var_three_name = 'SW_out'
            var_three_color = 'blue'
            var_four_name = 'LW_out'
            var_four_color = 'navy'
            units = 'w/m2'
            title = usage + ' SW and LW Radiation'

        elif code == 3:  # Potential vs Actual Incoming Shortwave Radiation
            var_one_name = 'SW_in'
            var_one_color = 'red'
            var_two_name = 'SW_pot_ASCE'
            var_two_color = 'black'
            var_three_name = 'null'
            var_three_color = 'black'
            var_four_name = 'null'
            var_four_color = 'black'
            units = 'w/m2'
            title = usage + ' Potential and Actual Incoming Shortwave Radiation'

        elif code == 4:  # Temperature
            var_one_name = 'Average Temperature'
            var_one_color = 'black'
            var_two_name = 'null'
            var_two_color = 'black'
            var_three_name = 'null'
            var_three_color = 'black'
            var_four_name = 'null'
            var_four_color = 'black'
            units = 'degrees Celsius'
            title = usage + var_one_name

        elif code == 5:  # Vapor Pressure and VPD
            var_one_name = 'Vapor Pressure'
            var_one_color = 'blue'
            var_two_name = 'Vapor Pressure Deficit'
            var_two_color = 'black'
            var_three_name = 'null'
            var_three_color = 'black'
            var_four_name = 'null'
            var_four_color = 'black'
            units = 'kPa'
            title = usage + var_one_name + ' and ' + var_two_name

        elif code == 6:  # Windspeed
            var_one_name = 'Windspeed'
            var_one_color = 'black'
            var_two_name = 'null'
            var_two_color = 'black'
            var_three_name = 'null'
            var_three_color = 'black'
            var_four_name = 'null'
            var_four_color = 'black'
            units = 'm/s'
            title = usage + var_one_name

        elif code == 7:  # Precipitation
            var_one_name = 'Precipitation'
            var_one_color = 'black'
            var_two_name = 'null'
            var_two_color = 'black'
            var_three_name = 'null'
            var_three_color = 'black'
            var_four_name = 'null'
            var_four_color = 'black'
            units = 'mm/day'
            title = usage + var_one_name

        elif code == 8:  # ET Comparison
            var_one_name = 'ET'
            var_one_color = 'black'
            var_two_name = 'ET_corr'
            var_two_color = 'blue'
            var_three_name = 'ET_user_corr'
            var_three_color = 'skyblue'
            var_four_name = 'null'
            var_four_color = 'black'
            units = 'mm/day'
            title = usage + ' Evapotraspiration'

        elif code == 81:  # ET Comparison no user_corr
            var_one_name = 'ET'
            var_one_color = 'black'
            var_two_name = 'ET_corr'
            var_two_color = 'blue'
            var_three_name = 'null'
            var_three_color = 'skyblue'
            var_four_name = 'null'
            var_four_color = 'black'
            units = 'mm/day'
            title = usage + ' Evapotraspiration'

        elif code == 9:  # LE comparison
            var_one_name = 'LE'
            var_one_color = 'black'
            var_two_name = 'LE_corr'
            var_two_color = 'blue'
            var_three_name = 'LE_user_corr'
            var_three_color = 'skyblue'
            var_four_name = 'null'
            var_four_color = 'black'
            units = 'w/m2'
            title = usage + ' Latent Heat Flux'

        elif code == 91:  # LE comparison no user_corr
            var_one_name = 'LE'
            var_one_color = 'black'
            var_two_name = 'LE_corr'
            var_two_color = 'blue'
            var_three_name = 'null'
            var_three_color = 'skyblue'
            var_four_name = 'null'
            var_four_color = 'black'
            units = 'w/m2'
            title = usage + ' Latent Heat Flux'

        elif code == 10:  # energy balance ratio comparison
            var_one_name = 'ebr'
            var_one_color = 'black'
            var_two_name = 'ebr_corr'
            var_two_color = 'blue'
            var_three_name = 'ebr_user_corr'
            var_three_color = 'skyblue'
            var_four_name = 'null'
            var_four_color = 'black'
            units = 'Ratio of Flux/Energy'
            title = usage + ' Energy Balance Ratio'

        elif code == 101:  # energy balance ratio comparison no user_corr
            var_one_name = 'ebr_raw'
            var_one_color = 'black'
            var_two_name = 'ebr_adj'
            var_two_color = 'blue'
            var_three_name = 'null'
            var_three_color = 'skyblue'
            var_four_name = 'null'
            var_four_color = 'black'
            units = 'Ratio of Flux/Energy'
            title = usage + ' Energy Balance Ratio'

        elif code == 11:  # energy balance scatter plot, raw vs corr
            var_one_name = 'Energy'
            var_one_color = 'null'
            var_two_name = 'Flux_raw'
            var_two_color = 'null'
            var_three_name = 'Flux_corr'
            var_three_color = 'null'
            var_four_name = 'null'
            var_four_color = 'null'
            units = 'null'
            title = usage + ' EBC, Raw vs. Corrected'

        elif code == 12:  # LE scatter plot, raw vs corr
            var_one_name = 'LE'
            var_one_color = 'null'
            var_two_name = 'LE_corr'
            var_two_color = 'null'
            var_three_name = 'null'
            var_three_color = 'null'
            var_four_name = 'null'
            var_four_color = 'null'
            units = 'null'
            title = usage + ' LE, Raw vs. Corrected'

        elif code == 13:  # ET scatter plot, raw vs corr
            var_one_name = 'ET'
            var_one_color = 'null'
            var_two_name = 'ET_corr'
            var_two_color = 'null'
            var_three_name = 'null'
            var_three_color = 'null'
            var_four_name = 'null'
            var_four_color = 'null'
            units = 'null'
            title = usage + ' ET, Raw vs. Corrected'

        else:
            raise ValueError('Unsupported code type {} passed to _generate_plot_features.'.format(code))

        if '%' in usage:
            units = '% difference'
        else:
            pass

        return units, title, var_one_name, var_one_color, var_two_name, var_two_color, var_three_name, \
            var_three_color, var_four_name, var_four_color

    def _generate_line_plot(self, x_size, y_size, dt_array, code, usage, var_one=None, var_two=None, var_three=None,
                            var_four=None, link_plot=None):
        """
            Creates a bokeh line plot for provided variables and links them if appropriate

            Parameters:
                x_size (int): x-axis size for plot
                y_size (int): y-axis size for plot
                dt_array (pandas dataframe): values for x-axis to label timestep, either daily or mean monthly
                code (int): code indicating what variables were passed
                usage (str): additional text indicating why plot is being created
                var_one (numpy array): 1D numpy array of first variable
                var_two (numpy array): 1D numpy array of second variable (if provided)
                var_three (numpy array): 1D numpy array of third variable (if provided)
                var_four (numpy array): 1D numpy array of fourth variable (if provided)
                link_plot (bokeh object): either nothing or the subplot we want to link x-axis with

            Returns:
                subplot (bokeh object): constructed figure to be assembled in aggregate plot
        """

        date_list = dt_array.tolist()
        source = ColumnDataSource(data=dict(date=date_list, v_one=var_one))
        empty_array = np.zeros(len(date_list))
        empty_array[:] = np.nan

        if var_two is None:
            source.add(empty_array, name='v_two')
        else:
            source.add(var_two, name='v_two')
        if var_three is None:
            source.add(empty_array, name='v_three')
        else:
            source.add(var_three, name='v_three')
        if var_four is None:
            source.add(empty_array, name='v_four')
        else:
            source.add(var_four, name='v_four')

        tooltips = [
            ('Index', '$index'),
            ('Date', '@date{%F}'),
            ('Value', '$y')]
        formatters = {'date': 'datetime'}

        (units, title, var_one_name, var_one_color, var_two_name, var_two_color,
         var_three_name, var_three_color, var_four_name, var_four_color) = self._generate_plot_features(code, usage)

        if usage == 'Monthly ':
            x_label = 'Monthly Timestep'
            x_axis_type = 'datetime'
        else:  # Anything else
            x_label = 'Daily Timestep'
            x_axis_type = 'datetime'

        if link_plot is None:  # No plot to link with
            subplot = figure(
                width=x_size, height=y_size, x_axis_type=x_axis_type,
                x_axis_label=x_label, y_axis_label=units, title=title,
                tools='pan, box_zoom, undo, reset, save')
        else:  # Plot is passed to link x-axis with
            subplot = figure(
                x_range=link_plot.x_range,
                width=x_size, height=y_size, x_axis_type=x_axis_type,
                x_axis_label=x_label, y_axis_label=units, title=title,
                tools='pan, box_zoom, undo, reset, save')

        subplot.line(x='date', y='v_one', line_color=var_one_color, legend=var_one_name, source=source)

        if var_two_name.lower() == 'null':
            pass
        else:
            subplot.line(x='date', y='v_two', line_color=var_two_color, legend=var_two_name, source=source)

        if var_three_name.lower() == 'null':
            pass
        else:
            subplot.line(x='date', y='v_three', line_color=var_three_color, legend=var_three_name, source=source)

        if var_four_name.lower() == 'null':
            pass
        else:
            subplot.line(x='date', y='v_four', line_color=var_four_color, legend=var_four_name, source=source)

        subplot.legend.location = 'bottom_left'
        subplot.sizing_mode = 'stretch_both'
        subplot.add_tools(HoverTool(tooltips=tooltips, formatters=formatters))

        return subplot

    def _generate_scatter_plot(self, x_size, y_size, x_label, y_label, code, usage,
                               var_one=None, var_two=None, var_three=None, link_plot=None):
        """
            Creates a bokeh scatter plot for provided variables and links them if appropriate. Var_one is the x axis,
            with var_two and var_three (if present) serving as y axis in comparison with one another

            Parameters:
                x_size (int): x-axis size for plot
                y_size (int): y-axis size for plot
                x_label (string): text to be displayed on bottom axis
                y_label (string): text to be displayed on vertical axis
                code (int): code indicating what variables were passed
                usage (str): additional text indicating why plot is being created
                var_one (numpy array): 1D numpy array of first variable
                var_two (numpy array): 1D numpy array of second variable (if provided)
                var_three (numpy array): 1D numpy array of third variable (if provided)
                link_plot (bokeh object): either nothing or the subplot we want to link x-axis with

            Returns:
                subplot (bokeh object): constructed scatterplot figure to be assembled in aggregate plot
        """
        data_length = len(var_one)

        (units, title, var_one_name, var_one_color, var_two_name, var_two_color,
         var_three_name, var_three_color, var_four_name, var_four_color) = self._generate_plot_features(code, usage)

        # Get max/mins from variables for the 1:1 line
        if var_one is not None:
            v1_max = np.nanmax(var_one)
            v1_min = np.nanmin(var_one)
        else:
            v1_max = np.nan
            v1_min = np.nan
        if var_two is not None:
            v2_max = np.nanmax(var_one)
            v2_min = np.nanmin(var_two)
        else:
            v2_max = np.nan
            v2_min = np.nan
        if var_three is not None:
            v3_max = np.nanmax(var_one)
            v3_min = np.nanmin(var_two)
        else:
            v3_max = np.nan
            v3_min = np.nan

        graph_max = np.nanmax([v1_max, v2_max, v3_max])
        graph_min = np.nanmin([v1_min, v2_min, v3_min])

        x = np.array([graph_min - 1, graph_max + 1])
        y = np.array([graph_min - 1, graph_max + 1])

        # have to loop elementwise to remove nans
        var_one_lstsq = []
        var_two_lstsq = []
        var_three_lstsq = []

        if var_three is not None:
            for i in range(0, data_length):
                if np.isnan(var_one[i]) or np.isnan(var_two[i]) or np.isnan(var_three[i]):
                    # do nothing, skipping over this entry because at least one component is a nan
                    pass
                else:
                    var_one_lstsq.append(var_one[i])
                    var_two_lstsq.append(var_two[i])
                    var_three_lstsq.append(var_three[i])

            var_one_lstsq = np.array(var_one_lstsq)
            var_two_lstsq = np.array(var_two_lstsq)
            var_three_lstsq = np.array(var_three_lstsq)

            # linalg.lstsq requires a column vector
            var_one_lstsq = var_one_lstsq[:, np.newaxis]

            slope_orig, _, _, _ = np.linalg.lstsq(
                var_one_lstsq, var_two_lstsq, rcond=None)
            slope_corr, _, _, _ = np.linalg.lstsq(
                var_one_lstsq, var_three_lstsq, rcond=None)

            if link_plot is None:  # No plot to link with
                subplot = figure(
                    width=x_size, height=y_size, x_axis_type='linear',
                    x_axis_label=x_label, y_axis_label=y_label, title=title,
                    tools='pan, box_zoom, undo, reset, hover, save')
            else:  # Plot is passed to link axis with
                subplot = figure(
                    x_range=link_plot.x_range, y_range=link_plot.y_range,
                    width=x_size, height=y_size, x_axis_type='linear',
                    x_axis_label=x_label, y_axis_label=y_label, title=title,
                    tools='pan, box_zoom, undo, reset, hover, save')

            subplot.line(x, y, line_color='black', legend='1:1 Line', line_dash='dashed')
            subplot.line(var_one, var_one * slope_orig, line_color='skyblue', legend='LSL_init, {:.3f}'
                         .format(float(slope_orig)))
            subplot.line(var_one, var_one * slope_corr, line_color='navy', legend='LSL_corr, {:.3f}'
                         .format(float(slope_corr)))

            subplot.triangle(var_one, var_two, size=10, color="lightsalmon", alpha=0.5, legend=var_two_name)
            subplot.circle(var_one, var_three, size=10, color="lightcoral", alpha=0.5, legend=var_three_name)
            subplot.legend.location = 'top_left'

        else:
            for i in range(0, data_length):
                if np.isnan(var_one[i]) or np.isnan(var_two[i]):
                    # do nothing, skipping over this entry because at least one component is a nan
                    pass
                else:
                    var_one_lstsq.append(var_one[i])
                    var_two_lstsq.append(var_two[i])

            var_one_lstsq = np.array(var_one_lstsq)
            var_two_lstsq = np.array(var_two_lstsq)

            # linalg.lstsq requires a column vector
            var_one_lstsq = var_one_lstsq[:, np.newaxis]

            slope_orig, _, _, _ = np.linalg.lstsq(
                var_one_lstsq, var_two_lstsq, rcond=None)

            if link_plot is None:  # No plot to link with
                subplot = figure(
                    width=x_size, height=y_size, x_axis_type='linear',
                    x_axis_label=x_label, y_axis_label=y_label, title=title,
                    tools='pan, box_zoom, undo, reset, hover, save')
            else:  # Plot is passed to link axis with
                subplot = figure(
                    x_range=link_plot.x_range, y_range=link_plot.y_range,
                    width=x_size, height=y_size, x_axis_type='linear',
                    x_axis_label=x_label, y_axis_label=y_label, title=title,
                    tools='pan, box_zoom, undo, reset, hover, save')

            subplot.line(x, y, line_color='black', legend='1:1 Line', line_dash='dashed')
            subplot.line(var_one, var_one * slope_orig, line_color='skyblue', legend='LSL, {:.3f}'
                         .format(float(slope_orig)))

            subplot.triangle(var_one, var_two, size=10, color="lightsalmon", alpha=0.5, legend=var_two_name)
            subplot.legend.location = 'top_left'

        subplot.sizing_mode = 'stretch_both'
        return subplot

    def _generate_scalable_plot(self, x_size, y_size, dt_array, units, title, usage, df, inv_map, link_plot=None):
        """
            Creates a bokeh line plot for as many of a certain type of variable are provided, such as soil moisture
            or ground heat flux sensors

            Parameters:
                x_size (int): x-axis size for plot
                y_size (int): y-axis size for plot
                dt_array (pandas dataframe): values for x-axis to label timestep, either daily or mean monthly
                usage (str): additional text indicating why plot is being created
                df (pandas dataframe): all of the variables that will be plotted on this graph
                inv_map (dictionary): map of
                link_plot (bokeh object): either nothing or the subplot we want to link x-axis with

            Returns:
                subplot (bokeh object): constructed figure to be assembled in aggregate plot
        """
        # Obtain original user names for display
        vars_to_plot = [col for col in df.columns]
        num_lines = len(vars_to_plot)
        color_list = Viridis256[0:-1:int(256/num_lines)]

        user_provided_names = []
        for el in vars_to_plot:
            user_provided_names.append(self.variables.get(el))

        if usage == 'Monthly ':
            x_label = 'Monthly Timestep'
            x_axis_type = 'datetime'
            title = usage + title
        else:  # Anything else
            x_label = 'Daily Timestep'
            x_axis_type = 'datetime'

        if link_plot is None:  # No plot to link with
            subplot = figure(
                width=x_size, height=y_size, x_axis_type=x_axis_type,
                x_axis_label=x_label, y_axis_label=units, title=title,
                tools='pan, box_zoom, undo, reset, hover, save')
        else:  # Plot is passed to link x-axis with
            subplot = figure(
                x_range=link_plot.x_range,
                width=x_size, height=y_size, x_axis_type=x_axis_type,
                x_axis_label=x_label, y_axis_label=units, title=title,
                tools='pan, box_zoom, undo, reset, hover, save')

        while len(vars_to_plot) > 0:
            var_name = vars_to_plot.pop(0)
            user_name = user_provided_names.pop(0)
            line_color = color_list.pop(0)

            subplot.line(dt_array, df[var_name], line_color=line_color, legend=user_name)
        subplot.legend.location = 'bottom_left'
        subplot.sizing_mode = 'stretch_both'
        return subplot

    def create_and_aggregate_plots(self, provided_vars, df, monthly_df):
        """
            Creates subplots for all of the different variable groupings that have been provided by the input data,
            determined by function _inventory_variables. Graphs will generate so long as one of the variables
            used is present

            Parameters:
                provided_vars (list): python list of all variables that have been provided
                df (pandas dataframe): dataframe of all variables, both provided and calculated
                monthly_df (pandas dataframe): dataframe of all variables for the monthly timestep

            Returns:
                fig (bokeh object): a grid of all of the constructed subplots
        """
        user_requested_columns = 1  # TODO: import number of columns from user
        x_size = 500
        y_size = 350
        plot_list = []
        monthly_plot_list = []

        # Plot surface balance components
        if ('Rn' in provided_vars) and ('G' in provided_vars) and \
                ('LE' in provided_vars) and ('H' in provided_vars):
            plot_surface_bal = self._generate_line_plot(
                    x_size, y_size, df.index, 1, '', df.Rn, df.LE, df.H, df.G)

            monthly_plot_surface_bal = self._generate_line_plot(
                    x_size, y_size, monthly_df.index, 1, 'Monthly ',
                    monthly_df.Rn, monthly_df.LE, monthly_df.H, monthly_df.G)
            plot_list.append(plot_surface_bal)
            monthly_plot_list.append(monthly_plot_surface_bal)
        else:
            print('\nSurface balance components graph missing all variables.')
            plot_surface_bal = None
            monthly_plot_surface_bal = None

        # Plot net radiation components
        if ('sw_in' in provided_vars) and ('sw_out' in provided_vars) and \
                ('lw_in' in provided_vars) and ('lw_out' in provided_vars):

            plot_net_rad = self._generate_line_plot(
                    x_size, y_size, df.index, 2, '', df.sw_in, df.lw_in, 
                    df.sw_out, df.lw_out)

            monthly_plot_net_rad = self._generate_line_plot(
                    x_size, y_size, monthly_df.index, 2, 'Monthly ',
                    monthly_df.sw_in, monthly_df.lw_in, monthly_df.sw_out, 
                    monthly_df.lw_out)

            plot_list.append(plot_net_rad)
            monthly_plot_list.append(monthly_plot_net_rad)
        else:
            print('\nNet Radiation components graph missing all variables.')
            plot_net_rad = None
            monthly_plot_net_rad = None

        # Plot temperature
        if 't_avg' in provided_vars:
            plot_temp = self._generate_line_plot(
                    x_size, y_size, df.index, 4, '', df.t_avg, None, None, None)
            monthly_plot_temp = self._generate_line_plot(
                    x_size, y_size, monthly_df.index, 4, 'Monthly ',
                    monthly_df.t_avg, None, None, None)
            plot_list.append(plot_temp)
            monthly_plot_list.append(monthly_plot_temp)
        else:
            print('\nTemperature graph missing all variables.')
            plot_temp = None
            monthly_plot_temp = None

        # Plot Vapor_Pres and VPD
        if ('vp' in provided_vars) and ('vpd' in provided_vars):
            plot_vapor_pres = self._generate_line_plot(
                    x_size, y_size, df.index, 5, '', df.vp, df.vpd, None, None)
            monthly_plot_vapor_pres = self._generate_line_plot(
                    x_size, y_size, monthly_df.index, 5, 'Monthly ', monthly_df.vp, monthly_df.vpd, None, None)
            plot_list.append(plot_vapor_pres)
            monthly_plot_list.append(monthly_plot_vapor_pres)

        else:
            print('\nVapor Pressure graph missing a variable.')
            plot_vapor_pres = None
            monthly_plot_vapor_pres = None

        # Plot windspeed
        if 'ws' in provided_vars:
            plot_windspeed = self._generate_line_plot(
                    x_size, y_size, df.index, 6, '', df.ws, None, None, None)
            monthly_plot_windspeed = self._generate_line_plot(
                    x_size, y_size, monthly_df.index, 6, 'Monthly ',
                    monthly_df.ws, None, None, None)
            plot_list.append(plot_windspeed)
            monthly_plot_list.append(monthly_plot_windspeed)
        else:
            print('\nWindspeed graph missing all variables.')
            plot_windspeed = None
            monthly_plot_windspeed = None

        # Plot precipitation
        if 'ppt' in provided_vars:
            plot_precip = self._generate_line_plot(
                    x_size, y_size, df.index, 7, '', df.ppt, None, None, None)
            monthly_plot_precip = self._generate_line_plot(
                    x_size, y_size, monthly_df.index, 7, 'Monthly ',
                    monthly_df.ppt, None, None, None)
            plot_list.append(plot_precip)
            monthly_plot_list.append(monthly_plot_precip)
        else:
            print('\nPrecipitation graph missing all variables.')
            plot_precip = None
            monthly_plot_precip = None

        # Plot ET
        if 'et' in provided_vars:
            # still plot if user did not have their own corrected data
            if not 'et_user_corr' in provided_vars:
                plot_et = self._generate_line_plot(
                    x_size, y_size, df.index, 81, '', df.et, df.et_corr, None)
                monthly_plot_et = self._generate_line_plot(x_size, y_size, 
                    monthly_df.index, 81, 'Monthly ', monthly_df.et,
                    monthly_df.et_corr, None)
            else:
                plot_et = self._generate_line_plot(
                        x_size, y_size, df.index, 8, '', df.et, df.et_corr,
                        df.et_user_corr, None)
                monthly_plot_et = self._generate_line_plot(x_size, y_size, 
                        monthly_df.index, 8, 'Monthly ', monthly_df.et,
                        monthly_df.et_corr, monthly_df.et_user_corr, None)
            plot_list.append(plot_et)
            monthly_plot_list.append(monthly_plot_et)
        else:
            print('\nEvapotranspiration graph missing a variable.')
            plot_et = None
            monthly_plot_et = None

        # Plot LE
        if 'LE' in provided_vars and 'LE_corr' in provided_vars:
            # still plot if no user corrected data
            if not 'LE_user_corr' in provided_vars:
                plot_le = self._generate_line_plot(
                    x_size, y_size, df.index, 91, '', df.LE, df.LE_corr, None)
                monthly_plot_le = self._generate_line_plot(
                    x_size, y_size, monthly_df.index, 91, 'Monthly ', 
                    monthly_df.LE, monthly_df.LE_corr, None)
            else:
                plot_le = self._generate_line_plot(
                    x_size, y_size, df.index, 9, '', df.LE, df.LE_corr, 
                    df.LE_user_corr, None)
                monthly_plot_le = self._generate_line_plot(
                    x_size, y_size, monthly_df.index, 9, 'Monthly ', 
                    monthly_df.LE, monthly_df.LE_corr, monthly_df.LE_user_corr, None)
            plot_list.append(plot_le)
            monthly_plot_list.append(monthly_plot_le)
        else:
            print('\nLE graph missing a variable.')
            plot_le = None
            monthly_plot_le = None

        # Plot energy balance ratios
        if ('ebr' in provided_vars) and ('ebr_corr' in provided_vars):
            # still plot if without user corrected data
            if not 'ebr_user_corr' in provided_vars:
                plot_ebr = self._generate_line_plot(
                    x_size, y_size, df.index, 101, '', df.ebr, df.ebr_corr,
                    None)
                monthly_plot_ebr = self._generate_line_plot(
                    x_size, y_size, monthly_df.index, 101, 'Monthly ', 
                    monthly_df.ebr, monthly_df.ebr_corr, None)

                # Create label of overall averages between EBR approaches
                avg_ebr_raw = (df.LE.mean() + df.H.mean()) / \
                        (df.Rn.mean() - df.G.mean())
                avg_ebr_corr = (df.LE_corr.mean() + df.H_corr.mean()) / \
                        (df.Rn.mean() - df.G.mean())

                ratio_averages_label = Label(x=70, y=225, x_units='screen', 
                        y_units='screen', 
                        text='EBR:{}\nEBR_corr:{}'\
                            .format(
                                str(round(avg_ebr_raw,3)), 
                                str(round(avg_ebr_corr,3))),
                        render_mode='css', border_line_color='black', 
                        border_line_alpha=0.25, background_fill_color='white', 
                        background_fill_alpha=1.0) 
                # create a copy of the label because bokeh doesnt want to assign 
                # the same label twice
                monthly_ratio_averages_label = Label(x=70, y=225, 
                        x_units='screen', y_units='screen', 
                        text='EBR:{}\nEBR_corr:{}'\
                            .format(
                                str(round(avg_ebr_raw,3)), 
                                str(round(avg_ebr_corr,3))),
                        render_mode='css', border_line_color='black', 
                        border_line_alpha=0.25, background_fill_color='white', 
                        background_fill_alpha=1.0)
            else:
                plot_ebr = self._generate_line_plot(
                    x_size, y_size, df.index, 10, '', df.ebr, df.ebr_corr,
                    df.ebr_user_corr, None)
                monthly_plot_ebr = self._generate_line_plot(
                    x_size, y_size, monthly_df.index, 10, 'Monthly ', 
                    monthly_df.ebr, monthly_df.ebr_corr,
                    monthly_df.ebr_user_corr, None)

                # Create label of overall averages between EBR approaches
                avg_ebr_raw = (df.LE.mean() + df.H.mean()) / \
                        (df.Rn.mean() - df.G.mean())
                avg_ebr_corr = (df.LE_corr.mean() + df.H_corr.mean()) / \
                        (df.Rn.mean() - df.G.mean())
                avg_ebr_user_corr = (df.LE_user_corr.mean() + df.H_user_corr.mean()) / \
                        (df.Rn.mean() - df.G.mean())

                ratio_averages_label = Label(x=70, y=225, x_units='screen', 
                        y_units='screen', 
                        text='EBR_raw:{:.3f}\nEBR_corr:{:.3f} \nEBR_user_corr:{:.3f}'\
                            .format(avg_ebr_raw, avg_ebr_corr, avg_ebr_user_corr),
                        render_mode='css', border_line_color='black', 
                        border_line_alpha=0.25, background_fill_color='white', 
                        background_fill_alpha=1.0) 
                # create a copy of the label because bokeh doesnt want to assign 
                # the same label twice
                monthly_ratio_averages_label = Label(x=70, y=225, 
                        x_units='screen', y_units='screen', 
                        text='EBR_raw:{:.3f}\nEBR_corr:{:.3f}\nEBR_user_corr:{:.3f}'\
                            .format(avg_ebr_raw, avg_ebr_corr, avg_ebr_user_corr),
                        render_mode='css', border_line_color='black', 
                        border_line_alpha=0.25, background_fill_color='white', 
                        background_fill_alpha=1.0)

            plot_ebr.add_layout(ratio_averages_label)
            monthly_plot_ebr.add_layout(monthly_ratio_averages_label)

            plot_list.append(plot_ebr)
            monthly_plot_list.append(monthly_plot_ebr)
        else:
            print('\nEnergy balance ratio graph missing a variable.')
            plot_ebr = None
            monthly_plot_ebr = None

        # Comparison plots of EBC
        if ('energy' in provided_vars) and ('flux' in provided_vars):
            plot_ebr_scatter = self._generate_scatter_plot(
                x_size, y_size, 'Energy (Net Radiation - G)', 'Flux (Latent Heat + Sensible Heat)', 11,
                '', var_one=df.energy, var_two=df.flux, var_three=df.flux_corr)

            monthly_plot_ebr_scatter = self._generate_scatter_plot(
                x_size, y_size, 'Energy (Net Radiation - G)', 'Flux (Latent Heat + Sensible Heat)', 11,
                'Monthly ', var_one=monthly_df.energy,
                var_two=monthly_df.flux, var_three=monthly_df.flux_corr)

            plot_list.append(plot_ebr_scatter)
            monthly_plot_list.append(monthly_plot_ebr_scatter)
        else:
            print('\nEnergy balance correlation graphs missing a variable.')
            plot_ebr_scatter = None
            monthly_plot_ebr_scatter = None

        # Correlation plots of LE vs LE_corr
        if ('LE' in provided_vars) and ('LE_corr' in provided_vars):
            plot_le_scatter = self._generate_scatter_plot(
                x_size, y_size, 'LE (w/m2)', 'LE_corr (w/m2)', 12,
                '', var_one=df.LE, var_two=df.LE_corr, var_three=None)

            monthly_plot_le_scatter = self._generate_scatter_plot(
                x_size, y_size, 'LE (w/m2)', 'LE_corr (w/m2)', 12,
                'Monthly ', var_one=monthly_df.LE, var_two=monthly_df.LE_corr, var_three=None)

            plot_list.append(plot_le_scatter)
            monthly_plot_list.append(monthly_plot_le_scatter)
        else:
            print('\n Latent Heat scatter plot missing a variable.')
            plot_le_scatter = None
            monthly_plot_le_scatter = None

        # Correlation plots of ET and ET_corr
        if ('et' in provided_vars) and ('et_corr' in provided_vars):
            plot_et_scatter = self._generate_scatter_plot(
                x_size, y_size, 'ET (mm)', 'ET_corr (mm)', 13,
                '', var_one=df.et, var_two=df.et_corr, var_three=None)

            monthly_plot_et_scatter = self._generate_scatter_plot(
                x_size, y_size, 'ET (mm)', 'ET_corr (mm)', 13,
                'Monthly ', var_one=monthly_df.et, var_two=monthly_df.et_corr, var_three=None)

            plot_list.append(plot_et_scatter)
            monthly_plot_list.append(monthly_plot_et_scatter)
        else:
            print('\n Latent Heat scatter plot missing a variable.')
            plot_et_scatter = None
            monthly_plot_et_scatter = None

        if 'theta_mean' in provided_vars:

            theta_cols = [col for col in df.columns if 'theta_' in col]
            theta_df = df[theta_cols].copy()
            monthly_theta_df = monthly_df[theta_cols].copy()
            multi_theta = self._generate_scalable_plot(
                x_size, y_size, df.index, 'cm', 'Soil Moisture', '', theta_df, self.inv_map, None)
            monthly_multi_theta = self._generate_scalable_plot(
                x_size, y_size, df.index, 'cm', 'Soil Moisture', 'Monthly ', monthly_theta_df, self.inv_map, None)

            plot_list.append(multi_theta)
            monthly_plot_list.append(monthly_multi_theta)
        else:
            print('\nSoil Moisture scalable plot missing a variable.')
            multi_theta = None
            monthly_multi_theta = None

        if 'g_mean' in provided_vars or 'G_mean' in provided_vars:

            g_cols = [col for col in df.columns if ((('g_' in col) or ('G' in col)) and not ('qc' in col))]
            g_df = df[g_cols].copy()
            monthly_g_df = monthly_df[g_cols].copy()
            multi_g = self._generate_scalable_plot(
                x_size, y_size, df.index,  'w/m_2', 'Soil Heat Flux', '', g_df, self.inv_map, None)
            monthly_multi_g = self._generate_scalable_plot(
                x_size, y_size, df.index, 'w/m_2', 'Soil Heat Flux', 'Monthly ', monthly_g_df, self.inv_map, None)

            plot_list.append(multi_g)
            monthly_plot_list.append(monthly_multi_g)
        else:
            print('\nSoil heat flux scalable plot missing a variable.')
            multi_g = None
            monthly_multi_g = None

        number_of_plots = len(plot_list)*2
        number_of_rows = ceil(number_of_plots / user_requested_columns)
        grid_of_plots = [[None] * user_requested_columns for i in range(number_of_rows)]

        for i in range(number_of_rows):
            for j in range(user_requested_columns):

                if len(plot_list) > 0:
                    grid_of_plots[i][j] = plot_list.pop(0)
                elif len(plot_list) == 0 and len(monthly_plot_list) > 0:
                    grid_of_plots[i][j] = monthly_plot_list.pop(0)
                else:
                    pass

        fig = gridplot(grid_of_plots, toolbar_location='left')


        return fig

    def generate_plots(self):
        """
        Create all of the graphs for provided variables, this is the function that actually runs all the code



        """
        if not self.out_dir.is_dir():
            self.out_dir.mkdir(parents=True, exist_ok=True)

        figure_path = self.out_dir.joinpath(
            '{}_plots.html'.format(self.site_id)
        )
        output_file(figure_path)

        # next lines rename columns to fluxdataqaqc names
        df = self._df.rename(columns=self.inv_map)
        monthly_df = self._monthly_df.rename(columns=self.inv_map)
        compound_fig = self.create_and_aggregate_plots(
                self.provided_vars, df, monthly_df)
        save(compound_fig)
        self.plots_created = True
        self.plot_file = figure_path

