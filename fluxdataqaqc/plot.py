# -*- coding: utf-8 -*-
"""
Plot fluxnet data using bokeh packages
"""

from .qaqc import QaQc
from bokeh.layouts import gridplot
from bokeh.models import Label
from bokeh.plotting import figure, output_file, save
from pathlib import Path
import numpy as np
from math import ceil


class Plot(object):
    """
    Takes a Data object from data.py, determines what variables are available, and makes bokeh plots for the different
    climate variables.

    """

    def __init__(self, qaqc=None):

        if isinstance(qaqc, QaQc):
            self._df = qaqc.df
            self._monthly_df = qaqc.monthly_df
            self.directory = qaqc.directory
            self.climate_file_name = qaqc.climate_file_name
            self.provided_vars = self._inventory_variables(self._df)
        elif qaqc is not None:
            raise TypeError("Must assign a fluxdataqaqc.qaqc.QaQc object")
        else:
            self._df = None

        self.plots_created = False

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
            title = usage + ' EBC, Raw vs. Fluxnet Corrected'

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

        data_length = len(var_one)

        (units, title, var_one_name, var_one_color, var_two_name, var_two_color,
         var_three_name, var_three_color, var_four_name, var_four_color) = self._generate_plot_features(code, usage)

        # Comparison plots of EBC
        x_label = 'Energy (Net Radiation - G)'
        y_label = 'Flux (Latent Heat + Sensible Heat)'

        # 1:1 line
        x = np.array([-80, 280])
        y = np.array([-80, 280])

        # have to loop elementwise to remove nans, TODO there is probably a better way to do this
        var_one_lstsq = []
        var_two_lstsq = []
        var_three_lstsq = []

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

        slope_orig, _, _, _ = np.linalg.lstsq(var_one_lstsq, var_two_lstsq)
        slope_corr, _, _, _ = np.linalg.lstsq(var_one_lstsq, var_three_lstsq)

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
        subplot.line(var_one, var_one*slope_orig, line_color='skyblue', legend='LSL_raw, {:.3f}'
                     .format(float(slope_orig)))
        subplot.line(var_one, var_one*slope_corr, line_color='navy', legend='LSL_corr, {:.3f}'
                     .format(float(slope_corr)))

        subplot.triangle(var_one, var_two, size=5, color="lightsalmon", alpha=0.5, legend=var_two_name)
        subplot.circle(var_one, var_three, size=5, color="lightcoral", alpha=0.5, legend=var_three_name)
        subplot.legend.location = 'top_left'
        return subplot

    def _generate_scalable_plot(self):
        pass

    def create_and_aggregate_plots(self, provided_vars, df, monthly_df):

        x_size = 500
        y_size = 350
        plot_list = []
        monthly_plot_list = []

        # Plot surface balance components
        if ('Rn' in provided_vars) and ('G' in provided_vars) and ('LE' in provided_vars) and ('H' in provided_vars):
            plot_surface_bal = self._generate_line_plot(x_size, y_size, df.index, 1, '', df.Rn, df.LE,
                                                    df.H, df.G)

            monthly_plot_surface_bal = self._generate_line_plot(x_size, y_size, monthly_df.index, 1, 'Monthly ',
                                                                monthly_df.Rn, monthly_df.LE,
                                                                monthly_df.H, monthly_df.G)
            plot_list.append(plot_surface_bal)
            monthly_plot_list.append(monthly_plot_surface_bal)
        else:
            print('\nSurface balance components graph missing a variable.')
            plot_surface_bal = None
            monthly_plot_surface_bal = None

        # Plot net radiation components
        if ('sw_in' in provided_vars) and ('sw_out' in provided_vars) and \
                ('lw_in' in provided_vars) and ('lw_out' in provided_vars):

            plot_net_rad = self._generate_line_plot(x_size, y_size, df.index, 2, '', df.sw_in, df.lw_in, df.sw_out,
                                                df.lw_out)

            monthly_plot_net_rad = self._generate_line_plot(x_size, y_size, monthly_df.index, 2, 'Monthly ',
                                                            monthly_df.sw_in, monthly_df.lw_in, monthly_df.sw_out,
                                                            monthly_df.lw_out)

            plot_list.append(plot_net_rad)
            monthly_plot_list.append(monthly_plot_net_rad)
        else:
            print('\nNet Radiation components graph missing a variable.')
            plot_net_rad = None
            monthly_plot_net_rad = None

        # Plot potential vs. measured inc_sw_rad
        if 'sw_in' in provided_vars:
            # TODO add rso to agg dict and include it here
            plot_sw_rad = self._generate_line_plot(x_size, y_size, df.index, 3, '', df.sw_in, df.rso, None)

            monthly_plot_sw_rad = self._generate_line_plot(x_size, y_size, monthly_df.index, 3, 'Monthly ',
                                                           monthly_df.sw_pot, monthly_df.sw_in, monthly_df.rso, None)
            plot_list.append(plot_sw_rad)
            monthly_plot_list.append(monthly_plot_sw_rad)
        else:
            print('\nMeasured vs Potential SW graph missing a variable.')
            plot_sw_rad = None
            monthly_plot_sw_rad = None

        # Plot temperature
        if 't_avg' in provided_vars:
            plot_temp = self._generate_line_plot(x_size, y_size, df.index, 4, '', df.t_avg, None, None, None)
            monthly_plot_temp = self._generate_line_plot(x_size, y_size, monthly_df.index, 4, 'Monthly ',
                                                         monthly_df.t_avg, None, None, None)
            plot_list.append(plot_temp)
            monthly_plot_list.append(monthly_plot_temp)
        else:
            print('\nTemperature graph missing a variable.')
            plot_temp = None
            monthly_plot_temp = None

        # Plot Vapor_Pres and VPD
        if ('vp' in provided_vars) and ('vpd' in provided_vars):
            plot_vapor_pres = self._generate_line_plot(x_size, y_size, df.index, 5, '', df.vp, df.vpd, None, None)
            monthly_plot_vapor_pres = self._generate_line_plot(x_size, y_size, monthly_df.index, 5, 'Monthly ',
                                                               monthly_df.vp, monthly_df.vpd, None, None)
            plot_list.append(plot_vapor_pres)
            monthly_plot_list.append(monthly_plot_vapor_pres)

        else:
            print('\nVapor Pressure graph missing a variable.')
            plot_vapor_pres = None
            monthly_plot_vapor_pres = None

        # Plot windspeed
        if 'ws' in provided_vars:
            plot_windspeed = self._generate_line_plot(x_size, y_size, df.index, 6, '', df.ws, None, None, None)
            monthly_plot_windspeed = self._generate_line_plot(x_size, y_size, monthly_df.index, 6, 'Monthly ',
                                                              monthly_df.ws, None, None, None)
            plot_list.append(plot_windspeed)
            monthly_plot_list.append(monthly_plot_windspeed)
        else:
            print('\nWindspeed graph missing a variable.')
            plot_windspeed = None
            monthly_plot_windspeed = None

        # Plot precipitation
        if 'ppt' in provided_vars:
            plot_precip = self._generate_line_plot(x_size, y_size, df.index, 7, '', df.ppt, None, None, None)
            monthly_plot_precip = self._generate_line_plot(x_size, y_size, monthly_df.index, 7, 'Monthly ',
                                                           monthly_df.ppt, None, None, None)
            plot_list.append(plot_precip)
            monthly_plot_list.append(monthly_plot_precip)
        else:
            print('\nPrecipitation graph missing a variable.')
            plot_precip = None
            monthly_plot_precip = None

        # Plot ET
        if 'et_reg' in provided_vars:
            plot_et = self._generate_line_plot(x_size, y_size, df.index, 8, '', df.et_reg, df.et_corr, df.et_adj, None)
            monthly_plot_et = self._generate_line_plot(x_size, y_size, monthly_df.index, 8, 'Monthly ',
                                                       monthly_df.et_reg, monthly_df.et_corr, monthly_df.et_adj, None)
            plot_list.append(plot_et)
            monthly_plot_list.append(monthly_plot_et)
        else:
            print('\nEvapotranspiration graph missing a variable.')
            plot_et = None
            monthly_plot_et = None

        # Plot LE
        if 'LE' in provided_vars:
            plot_le = self._generate_line_plot(x_size, y_size, df.index, 9, '', df.LE, df.LE_corr, df.LE_adj, None)
            monthly_plot_le = self._generate_line_plot(x_size, y_size, monthly_df.index, 9, 'Monthly ', monthly_df.LE,
                                                       monthly_df.LE_corr, monthly_df.LE_adj, None)
            plot_list.append(plot_le)
            monthly_plot_list.append(monthly_plot_le)
        else:
            print('\nLE graph missing a variable.')
            plot_le = None
            monthly_plot_le = None

        # Plot energy balance ratios
        if ('Rn' in provided_vars) and ('G' in provided_vars) and ('LE' in provided_vars) and ('H' in provided_vars):
            # TODO redundant condition, organize better
            plot_ebr = self._generate_line_plot(x_size, y_size, df.index, 10, '', df.ebc_reg, df.ebc_corr, df.ebc_adj,
                                                None)
            monthly_plot_ebr = self._generate_line_plot(x_size, y_size, monthly_df.index, 10, 'Monthly ',
                                                        monthly_df.ebc_reg, monthly_df.ebc_corr, monthly_df.ebc_adj,
                                                        None)

            # Create label of overall averages between EBR approaches
            avg_ebr_raw = (df.LE.mean() + df.H.mean()) / (df.Rn.mean() - df.G.mean())
            avg_ebr_corr = (df.LE_corr.mean() + df.H_corr.mean()) / (df.Rn.mean() - df.G.mean())
            avg_ebr_adj = (df.LE_adj.mean() + df.H_adj.mean()) / (df.Rn.mean() - df.G.mean())

            ratio_averages_label = Label(x=70, y=225, x_units='screen', y_units='screen',
                                         text='EBR_raw: {:.3f} \nEBR_corr: {:.3f} \n EBR_adj: {:.3f}'
                                         .format(avg_ebr_raw, avg_ebr_corr, avg_ebr_adj), render_mode='css',
                                         border_line_color='black', border_line_alpha=0.25,
                                         background_fill_color='white', background_fill_alpha=1.0)

            # create a copy of the label because bokeh doesnt want to assign the same label twice
            monthly_ratio_averages_label = Label(x=70, y=225, x_units='screen', y_units='screen',
                                                 text='EBR_raw: {:.3f} \nEBR_corr: {:.3f} \n EBR_adj: {:.3f}'
                                                 .format(avg_ebr_raw, avg_ebr_corr, avg_ebr_adj), render_mode='css',
                                                 border_line_color='black', border_line_alpha=0.25,
                                                 background_fill_color='white', background_fill_alpha=1.0)

            plot_ebr.add_layout(ratio_averages_label)
            monthly_plot_ebr.add_layout(monthly_ratio_averages_label)

            plot_list.append(plot_ebr)
            monthly_plot_list.append(monthly_plot_ebr)
        else:
            print('\nEnergy balance ratio graph missing a variable.')
            plot_ebr = None
            monthly_plot_ebr = None

        if ('energy' in provided_vars) and ('flux' in provided_vars):
            plot_ebc_corr = self._generate_scatter_plot(x_size, y_size, 11, '', var_one=df.energy, var_two=df.flux,
                                                        var_three=df.flux_corr)
            plot_ebc_adj = self._generate_scatter_plot(x_size, y_size, 12, '', var_one=df.energy, var_two=df.flux,
                                                       var_three=df.flux_adj)

            monthly_plot_ebc_corr = self._generate_scatter_plot(x_size, y_size, 11, 'Monthly ',
                                                                var_one=monthly_df.energy, var_two=monthly_df.flux,
                                                                var_three=monthly_df.flux_corr)
            monthly_plot_ebc_adj = self._generate_scatter_plot(x_size, y_size, 12, 'Monthly ',
                                                               var_one=monthly_df.energy, var_two=monthly_df.flux,
                                                               var_three=monthly_df.flux_adj)
            plot_list.append(plot_ebc_corr)
            monthly_plot_list.append(monthly_plot_ebc_corr)
            plot_list.append(plot_ebc_adj)
            monthly_plot_list.append(monthly_plot_ebc_adj)
        else:
            print('\nEnergy balance correlation graphs missing a variable.')
            plot_ebc_corr = None
            plot_ebc_adj = None
            monthly_plot_ebc_corr = None
            monthly_plot_ebc_adj = None

        number_of_plots = len(plot_list)*2
        number_of_rows = ceil(number_of_plots / 3)
        grid_of_plots = [[None] * 3 for i in range(number_of_rows)]

        for i in range(number_of_rows):
            for j in range(3):

                if len(plot_list) > 0:
                    grid_of_plots[i][j] = plot_list.pop(0)
                elif len(plot_list) == 0 and len(monthly_plot_list) > 0:
                    grid_of_plots[i][j] = monthly_plot_list.pop(0)
                else:
                    pass

        fig = gridplot(grid_of_plots, toolbar_location='left')
        # fig = gridplot([[plot_surface_bal, plot_net_rad, plot_sw_rad],
        #                 [plot_temp, plot_windspeed, plot_precip],
        #                 [plot_et, plot_le],
        #                 [plot_ebr, plot_ebc_corr, plot_ebc_adj],
        #                 [monthly_plot_surface_bal, monthly_plot_net_rad, monthly_plot_sw_rad],
        #                 [monthly_plot_temp, monthly_plot_windspeed, monthly_plot_precip],
        #                 [monthly_plot_et, monthly_plot_le],
        #                 [monthly_plot_ebr, monthly_plot_ebc_corr, monthly_plot_ebc_adj]], toolbar_location="left")

        return fig

    def generate_plots(self):
        """
        Create all of the graphs for provided variables,

        """
        figure_path = self.directory.joinpath(self.climate_file_name + '_figure.html')
        output_file(figure_path)

        compound_fig = self.create_and_aggregate_plots(self.provided_vars, self._df, self._monthly_df)
        save(compound_fig)
        self.plots_created = True

