# -*- coding: utf-8 -*-
"""
Plot fluxnet data using bokeh packages
"""

from .qaqc import QaQc
from bokeh.layouts import gridplot
from bokeh.plotting import figure, output_file, save
import numpy as np


class Plot(object):
    """
    Takes a Data object from data.py, determines what variables are available, and makes bokeh plots for the different
    climate variables.

    """

    def __init__(self, qaqc=None):

        if isinstance(qaqc, QaQc):
            self._df = qaqc.df
            self._monthly_df = qaqc.monthly_df
        elif qaqc is not None:
            raise TypeError("Must assign a fluxdataqaqc.qaqc.QaQc object")
        else:
            self._df = None

        self.provided_vars = self._inventory_variables(self._df)
        self.aggregate_plot = self._create_plots(self._df, self._monthly_df)

        output_file('example_output.html')
        save(self.aggregate_plot)

    def _inventory_variables(self, data):
        """
        Checks to see what variables are present in dataframe, here defined as if there are any non-nan values.
        A variable is given a flag of 0 if it IS missing, a flag of 1 indicates that it IS NOT missing

        :param data: pandas dataframe
        :return: A dictionary of flags for what variables are present
        """

        self.provided_vars = {}
        for variable in data.columns:
            observations = np.array(data[variable])
            if np.isnan(observations).all():
                self.provided_vars.update({variable: 0})
            else:
                self.provided_vars.update({variable: 1})

        return self.provided_vars

    def _generate_plot_features(self, code, usage=''):
        """
            Generates plot features depending on what code is passed

            Parameters:
                code : integer code passed by main script that indicates what type of data has been passed
                usage : string indicating why this plot is being created, it may be blank.

            Returns:
                units : string of units for passed variable
                title : string for title of plot
                var_one_name : string of first variable name
                var_one_color : string of color code to use for plotting variable one
                var_two_name : string of second variable name
                var_two_color : string of color code to use for plotting variable two
                var_three_name : string of third variable name
                var_three_color : string of color code to use for plotting variable three
                var_four_name : string of fourth variable name
                var_four_color : string of color code to use for plotting variable four
        """

        if code == 1:  # First graph, four variables
            var_one_name = 'Net_Rad'
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
            var_one_name = 'Incoming SW'
            var_one_color = 'red'
            var_two_name = 'Incoming LW'
            var_two_color = 'darkred'
            var_three_name = 'Outgoing SW'
            var_three_color = 'blue'
            var_four_name = 'Outgoing LW'
            var_four_color = 'navy'
            units = 'w/m2'
            title = usage + ' SW and LW Radiation'

        elif code == 3:  # Potential vs Actual Incoming Shortwave Radiation
            var_one_name = 'Incoming Shortwave Potential'
            var_one_color = 'black'
            var_two_name = 'Incoming Shortwave Actual'
            var_two_color = 'red'
            var_three_name = 'ASCE Rso_a'
            var_three_color = 'blue'
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
            var_one_name = 'ET_raw'
            var_one_color = 'black'
            var_two_name = 'ET_corr'
            var_two_color = 'blue'
            var_three_name = 'ET_adj'
            var_three_color = 'skyblue'
            var_four_name = 'null'
            var_four_color = 'black'
            units = 'mm/day'
            title = usage + ' Daily Evapotraspiration'

        elif code == 9:  # LE comparison
            var_one_name = 'LE_raw'
            var_one_color = 'black'
            var_two_name = 'LE_corr'
            var_two_color = 'blue'
            var_three_name = 'LE_adj'
            var_three_color = 'skyblue'
            var_four_name = 'null'
            var_four_color = 'black'
            units = 'w/m2'
            title = usage + ' Latent Heat Flux'

        elif code == 10:  # energy balance ratio comparison
            var_one_name = 'ebr_raw'
            var_one_color = 'black'
            var_two_name = 'ebr_corr'
            var_two_color = 'blue'
            var_three_name = 'ebr_adj'
            var_three_color = 'skyblue'
            var_four_name = 'null'
            var_four_color = 'black'
            units = 'Ratio of Flux/Energy'
            title = usage + ' Energy Balance Ratio'

        elif code == 11:  # energy balance scatter plot, raw vs fluxnet Corr
            var_one_name = 'Energy'
            var_one_color = 'null'
            var_two_name = 'Flux_raw'
            var_two_color = 'null'
            var_three_name = 'Flux_corr'
            var_three_color = 'null'
            var_four_name = 'null'
            var_four_color = 'null'
            units = 'null'
            title = usage + ' EBC, Raw vs. Fluxnet'

        elif code == 12:  # energy balance scatter plot, raw vs Bowen Corr
            var_one_name = 'Energy'
            var_one_color = 'null'
            var_two_name = 'Flux_raw'
            var_two_color = 'null'
            var_three_name = 'Flux_adj'
            var_three_color = 'null'
            var_four_name = 'null'
            var_four_color = 'null'
            units = 'null'
            title = usage + ' EBC, Raw vs. Bowen Adjusted'

        else:
            raise ValueError('Unsupported code type {} passed to generate_plot_features.'.format(code))

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
                x_size : x-axis size for plot
                y_size : y-axis size for plot
                dt_array : values for x-axis to label timestep, either daily or mean monthly
                code : integer indicating what variables were passed
                usage : additional string indicating why plot is being created
                var_one : 1D numpy array of first variable
                var_two : 1D numpy array of second variable
                var_three : 1D numpy array of third variable
                var_four : 1D numpy array of fourth variable
                link_plot : either nothing or the plot we want to link x-axis with

            Returns:
                subplot : constructed figure
        """
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
                tools='pan, box_zoom, undo, reset, hover, save')
        else:  # Plot is passed to link x-axis with
            subplot = figure(
                x_range=link_plot.x_range,
                width=x_size, height=y_size, x_axis_type=x_axis_type,
                x_axis_label=x_label, y_axis_label=units, title=title,
                tools='pan, box_zoom, undo, reset, hover, save')

        subplot.line(dt_array, var_one, line_color=var_one_color, legend=var_one_name)

        if var_two_name.lower() == 'null':
            pass
        else:
            subplot.line(dt_array, var_two, line_color=var_two_color, legend=var_two_name)

        if var_three_name.lower() == 'null':
            pass
        else:
            subplot.line(dt_array, var_three, line_color=var_three_color, legend=var_three_name)

        if var_four_name.lower() == 'null':
            pass
        else:
            subplot.line(dt_array, var_four, line_color=var_four_color, legend=var_four_name)

        subplot.legend.location = 'bottom_left'

        return subplot

    def _generate_scatter_plot(self, x_size, y_size, code, usage, var_one=None, var_two=None, var_three=None,
                               link_plot=None):

        (units, title, var_one_name, var_one_color, var_two_name, var_two_color,
         var_three_name, var_three_color, var_four_name, var_four_color) = self._generate_plot_features(code, usage)

        # Comparison plots of EBC
        x_label = 'Energy (Net Radiation - G)'
        y_label = 'Flux (Latent Heat - Sensible Heat)'

        # 1:1 line
        x = np.array([-80, 280])
        y = np.array([-80, 280])

        # linalg.lstsq requires a column vector
        # TODO test this more
        var_one_lstsq = var_one[:, np.newaxis]
        #slope_orig, _, _, _ = np.linalg.lstsq(var_one_lstsq, var_two)
        #slope_corr, _, _, _ = np.linalg.lstsq(var_one_lstsq, var_three)

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

        subplot.line(x, y, line_color='black', legend='1:1 Line')
        #subplot.line(var_one, var_one*slope_orig, line_color='orange', legend='LS Line, Raw, {}'.format(slope_orig))
        #subplot.line(var_one, var_one*slope_corr, line_color='navy', legend='LS Line, Corr, {}'.format(slope_corr))

        subplot.triangle(var_one, var_two, size=5, color="lightsalmon", alpha=0.5, legend=var_two_name)
        subplot.circle(var_one, var_three, size=5, color="lightcoral", alpha=0.5, legend=var_three_name)
        subplot.legend.location = 'bottom_left'
        return subplot

    def _create_plots(self, df, monthly_df):

        # TODO put if statements handling missing variables,

        x_size = 500
        y_size = 350

        # Create all the daily plots

        # Plot surface balance components
        plot_surface_bal = self._generate_line_plot(x_size, y_size, df.index, 1, '', df.Rn, df.LE,
                                                    df.H, df.G)
        # Plot net radiation components
        plot_net_rad = self._generate_line_plot(x_size, y_size, df.index, 2, '', df.sw_in, df.lw_in, df.sw_out,
                                                df.lw_out)
        # Plot potential vs. measured inc_sw_rad
        plot_sw_rad = self._generate_line_plot(x_size, y_size, df.index, 3, '', df.sw_pot, df.sw_in, df.rso, None)
        # Plot temperature
        plot_temp = self._generate_line_plot(x_size, y_size, df.index, 4, '', df.t_avg, None, None, None)
        # Plot Vapor_Pres and VPD
        # plot_vapor_pres = self._generate_line_plot(x_size, y_size, df.index, 5, '', df.vp, df.vpd, None, None)
        # Plot windspeed
        plot_windspeed = self._generate_line_plot(x_size, y_size, df.index, 6, '', df.ws, None, None, None)
        # Plot precipitation
        plot_precip = self._generate_line_plot(x_size, y_size, df.index, 7, '', df.ppt, None, None, None)
        # Plot ET
        plot_et = self._generate_line_plot(x_size, y_size, df.index, 8, '', df.et_reg, df.et_corr, df.et_adj, None)
        # Plot le
        plot_le = self._generate_line_plot(x_size, y_size, df.index, 9, '', df.LE, df.LE_corr, df.LE_adj,
                                           None)
        # Plot energy balance ratios
        plot_ebr = self._generate_line_plot(x_size, y_size, df.index, 10, '', df.ebc_reg, df.ebc_corr, df.ebc_adj,
                                            None)
        plot_ebc_corr = self._generate_scatter_plot(x_size, y_size, 11, '', var_one=df.energy, var_two=df.flux,
                                                    var_three=df.flux_corr)
        plot_ebc_adj = self._generate_scatter_plot(x_size, y_size, 12, '', var_one=df.energy, var_two=df.flux,
                                                   var_three=df.flux_adj)

        # Create all the monthly plots
        # Plot surface balance components
        monthly_plot_surface_bal = self._generate_line_plot(x_size, y_size, monthly_df.index, 1, 'Monthly ',
                                                            monthly_df.Rn, monthly_df.LE,
                                                            monthly_df.H, monthly_df.G)
        # Plot net radiation components
        monthly_plot_net_rad = self._generate_line_plot(x_size, y_size, monthly_df.index, 2, 'Monthly ',
                                                        monthly_df.sw_in, monthly_df.lw_in, monthly_df.sw_out,
                                                        monthly_df.lw_out)
        # Plot potential vs. measured inc_sw_rad
        monthly_plot_sw_rad = self._generate_line_plot(x_size, y_size, monthly_df.index, 3, 'Monthly ',
                                                       monthly_df.sw_pot, monthly_df.sw_in, monthly_df.rso, None)
        # Plot temperature
        monthly_plot_temp = self._generate_line_plot(x_size, y_size, monthly_df.index, 4, 'Monthly ', monthly_df.t_avg,
                                                     None, None, None)
        # Plot Vapor_Pres and VPD
        # monthly_plot_vapor_pres = self._generate_line_plot(x_size, y_size, monthly_df.index, 5, 'Monthly ',
        #                                                   monthly_df.vp, monthly_df.vpd, None, None)
        # Plot windspeed
        monthly_plot_windspeed = self._generate_line_plot(x_size, y_size, monthly_df.index, 6, 'Monthly ',
                                                          monthly_df.ws, None, None, None)
        # Plot precipitation
        monthly_plot_precip = self._generate_line_plot(x_size, y_size, monthly_df.index, 7, 'Monthly ', monthly_df.ppt,
                                                       None, None, None)
        # Plot ET
        monthly_plot_et = self._generate_line_plot(x_size, y_size, monthly_df.index, 8, 'Monthly ', monthly_df.et_reg,
                                                   monthly_df.et_corr, monthly_df.et_adj, None)
        # Plot le
        monthly_plot_le = self._generate_line_plot(x_size, y_size, monthly_df.index, 9, 'Monthly ', monthly_df.LE,
                                                   monthly_df.LE_corr, monthly_df.LE_adj,
                                                   None)
        # Plot energy balance ratios
        monthly_plot_ebr = self._generate_line_plot(x_size, y_size, monthly_df.index, 10, 'Monthly ',
                                                    monthly_df.ebc_reg, monthly_df.ebc_corr, monthly_df.ebc_adj,
                                                    None)
        monthly_plot_ebc_corr = self._generate_scatter_plot(x_size, y_size, 11, 'Monthly ', var_one=monthly_df.energy,
                                                            var_two=monthly_df.flux,
                                                            var_three=monthly_df.flux_corr)
        monthly_plot_ebc_adj = self._generate_scatter_plot(x_size, y_size, 12, 'Monthly ', var_one=monthly_df.energy,
                                                           var_two=monthly_df.flux,
                                                           var_three=monthly_df.flux_adj)

        fig = gridplot([[plot_surface_bal, plot_net_rad, plot_sw_rad],
                        [plot_temp, plot_windspeed, plot_precip],
                        [plot_et, plot_le],
                        [plot_ebr, plot_ebc_corr, plot_ebc_adj],
                        [monthly_plot_surface_bal, monthly_plot_net_rad, monthly_plot_sw_rad],
                        [monthly_plot_temp, monthly_plot_windspeed, monthly_plot_precip],
                        [monthly_plot_et, monthly_plot_le],
                        [monthly_plot_ebr, monthly_plot_ebc_corr, monthly_plot_ebc_adj]], toolbar_location="left")

        return fig





