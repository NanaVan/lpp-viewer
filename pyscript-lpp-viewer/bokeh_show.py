#!/usr/bin/env python3
# -*- coding:utf-8 -*-

from bokeh.plotting import figure, curdoc, show
from bokeh.models import ColumnDataSource, DataTable, TableColumn, ColorBar, LogColorMapper, NumericInput, AutocompleteInput, Button, Div, HoverTool, InlineStyleSheet, Checkbox, TabPanel, Tabs, LabelSet
from bokeh.events import ButtonClick
from bokeh.layouts import layout, row, column
from bokeh.palettes import YlOrRd9, YlGnBu9

import numpy as np
from iid import IID

class Bokeh_show():
    '''
    A bokeh to show the yield/weight heatmap with corresponding revolution and spectrum figure
    '''
    def __init__(self, lppion, cen_freq, span, gamma_t, delta_Brho_over_Brho, gamma_setting, delta_v_over_v=1e-6, L_CSRe=128.8):
        '''
        extract all the secondary fragments and their respective yields calculated by LISE++
        (including Mass, Half-life, Yield of all the fragments)
        lppion:     LISE++ output file to be loaded
        cen_freq:   center frequency of the spectrum in MHz
        span:       span of the spectrum in kHz
        n_peak:     number of peaks to be identified
        L_CSRe:     circumference of CSRe in m, default value 128.8
        '''
        self.iid = IID(lppion, cen_freq, span, gamma_t, delta_Brho_over_Brho, gamma_setting, delta_v_over_v, L_CSRe, False)
        self.panel_control()
        self._initial()
        self._initial_RevTime_spectrum()
        self.iid.calc_ecooler_peak()
        self._initial_ec_on()
        self.tabs_main = Tabs(tabs=[self.tabpanel_ec_off, self.tabpanel_ec_on, self.tabpanel_TOF])
        self.tabs_yield = Tabs(tabs=[self.tabpanel_yield_ec_off, self.tabpanel_yield_ec_on, self.tabpanel_yield_TOF])
        self.p_spectrum_default_linear.visible = False
        self.p_spectrum_cooler_linear.visible = False
        self.p_spectrum_TOF_linear.visible = False
        self.labels_default.visible = False
        self.labels_cooler.visible = False
        self.labels_TOF.visible = False
        print('Bokeh: initial start')

    def _log(self, status):
        '''
        send the status of the process
        '''
        self.div_log.text = '{:}'.format(status)

    def _show(self):
        '''
        return Bokeh_dislay
        
        yield heatmap:
        total yield of the ions corresponding to the LISE++ file, display in nuclei map
        x: N,
        y: Z,
        z: total yield, including bare, H-like, He-like, etc. but not including isomers
        tips:
            element, Z, N, isomer numbers, yields of bare, H-like, He-like, etc.


        simulation spectrum and table
        simulation spectrum of the ions corresponding to the LISE++ file and setting
        x: frequency [kHz],
        y: weight 
        tips:
            ion(isometric_state), peak location, weight, yield, harmonic, revolution frequency, half life
        
        '''
        print('Bokeh: initial complete')
        self._log('Bokeh: initial complete')
        return column([row([column([row(self.input_cen_freq, self.input_span), row([self.input_L_CSRe, self.input_delta_Brho_over_Brho]), row([self.input_gamma_t, self.input_alpha_p]), self.div_log, row([ self.input_gamma_setting, self.input_m_over_q]), row([self.input_delta_v_over_v, self.button_set]), self.checkbox_ec_on, row([self.input_Brho, self.input_peakloc, self.button_calibrate]), self.checkbox_Brho_input, row([self.input_ion, self.button_find_ion, self.button_reset_ion]), self.checkbox_TOF_on, row([self.input_show_threshold, self.checkbox_log_or_linear, self.checkbox_labels_on])]), self.tabs_yield]), self.tabs_main])

    def _wrap_data(self, data_type):
        '''
        data_type: 1 for ISOCHRONOUSION, 0 for ECOOLERION, -1 for Null, -2 for TOFION
        return data for plot and table
        '''
        if data_type > -1:
            if data_type:
                result = self.iid.cur.execute("SELECT ION, ELEMENT, N, Z, ISOMERIC, PEAKLOC, PEAKSIG, HARMONIC, REVFREQ, HALFLIFE, YIELD, TOTALYIELD, WEIGHT, PEAKMAX, GAMMA FROM ISOCHRONOUSION WHERE PEAKMAX>=?", (self.input_show_threshold.value,)).fetchall()
            else:
                result = self.iid.cur.execute("SELECT ION, ELEMENT, N, Z, ISOMERIC, PEAKLOC, PEAKSIG, HARMONIC, REVFREQ, HALFLIFE, YIELD, TOTALYIELD, WEIGHT, PEAKMAX, PSEUDOGAMMA FROM ECOOLERION WHERE PEAKMAX>=?", (self.input_show_threshold.value,)).fetchall()
            data = {
                'ion': [item[0] for item in result],
                'element': [item[1] for item in result],
                'N': [item[2] for item in result],
                'Z': [item[3] for item in result],
                'isometric': [item[4] for item in result],
                'peak_loc': [item[5] for item in result],
                'peak_sig': [item[6] for item in result],
                'harmonic': [item[7] for item in result],
                'rev_freq': [item[8] for item in result],
                'half_life': [item[9] for item in result],
                'yield': [item[10] for item in result],
                'total_yield': [item[11] for item in result],
                'weight': [item[12] for item in result],
                'peak_max': [item[13] for item in result],
                'ion_label': [item[0]+'('+item[4]+')' for item in result]
                }
            if len(data['ion']) == 0:
                data['xs'] = []
                data['ys'] = []
                data['color'] = []
            else:
                data['xs'] = [np.concatenate([peak_loc-10**np.linspace(2,-8,21), np.array([peak_loc]), peak_loc+10**np.linspace(-8,2,21)]) for peak_loc in data['peak_loc']]
                data['ys'] = [amp / sig / np.sqrt(2*np.pi) * np.exp(-np.concatenate([-10**np.linspace(2,-8,21), np.array([0]), 10**np.linspace(-8,2,21)])**2/ 2 / sig**2) for amp, sig in zip(data['weight'], data['peak_sig'])]
                yield_top = int(np.log10(np.max(data['total_yield'])))
                if data_type:
                    data['color'] = [YlOrRd9[yield_top-int(np.log10(item))] if yield_top-int(np.log10(item))<=8 else YlOrRd9[-1] for item in data['total_yield']]
                else:
                    data['color'] = [YlGnBu9[yield_top-int(np.log10(item))] if yield_top-int(np.log10(item))<=8 else YlGnBu9[-1] for item in data['total_yield']]
            if data_type:
                data['gamma'] = [item[14] for item in result]
            else:
                data['pseudo_gamma'] = [item[14] for item in result]
        elif data_type < -1:
            result = self.iid.cur.execute("SELECT ION, ELEMENT, N, Z, ISOMERIC, PEAKSIG, REVTIME, HALFLIFE, YIELD, TOTALYIELD, GAMMA, PEAKMAX FROM TOFION WHERE PEAKMAX>=?", (self.input_show_threshold.value,)).fetchall()
            data = {
                'ion': [item[0] for item in result],
                'element': [item[1] for item in result],
                'N': [item[2] for item in result],
                'Z': [item[3] for item in result],
                'isometric': [item[4] for item in result],
                'peak_sig': [item[5] for item in result],
                'rev_time': [item[6] for item in result],
                'half_life': [item[7] for item in result],
                'yield': [item[8] for item in result],
                'total_yield': [item[9] for item in result],
                'gamma': [item[10] for item in result],
                'peak_max': [item[11] for item in result],
                'ion_label': [item[0]+'('+item[4]+')' for item in result]
                }
            if len(data['ion']) == 0:
                data['xs'] = []
                data['ys'] = []
            else:
                data['xs'] = [np.concatenate([rev_time-10**np.linspace(1,-8,31), np.array([rev_time]), rev_time+10**np.linspace(-8,1,31)]) for rev_time in data['rev_time']]
                data['ys'] = [amp / sig / np.sqrt(2*np.pi) * np.exp(-np.concatenate([-10**np.linspace(1,-8,31), np.array([0]), 10**np.linspace(-8,1,31)])**2/ 2 / sig**2) for amp, sig in zip(data['yield'], data['peak_sig'])]
                yield_top = int(np.log10(np.max(data['total_yield'])))
                if data_type:
                    data['color'] = [YlOrRd9[yield_top-int(np.log10(item))] if yield_top-int(np.log10(item))<=8 else YlOrRd9[-1] for item in data['total_yield']]
                else:
                    data['color'] = [YlGnBu9[yield_top-int(np.log10(item))] if yield_top-int(np.log10(item))<=8 else YlGnBu9[-1] for item in data['total_yield']]
        else:
            data = {'xs':[], 'ys':[], 'ion': [], 'element':[], 'N':[], 'Z':[], 'isometric': [], 'peak_loc': [], 'peak_sig': [], 'harmonic': [], 'rev_freq': [], 'rev_time': [], 'half_life': [], 'yield': [], 'total_yield': [], 'weight': [], 'gamma': [], 'pseudo_gamma': []}
        return data


    def _initial(self):
        def selected_ion(attr, old, new):
            try:
                self.temp_ion = self.spectrum_source.data['ion'][new[0]]
                self.temp_isometric_state = self.spectrum_source.data['isometric'][new[0]]
                self.temp_harmonic = self.spectrum_source.data['harmonic'][new[0]]
                self.temp_gamma = self.spectrum_source.data['gamma'][new[0]]
                self.input_gamma_setting.value = self.temp_gamma
                print("{:}({:}), γ: {:.3f}".format(self.temp_ion, self.temp_isometric_state, self.temp_gamma))
                self._log("Ion Selected: {:}({:}), γ: {:.3f}".format(self.temp_ion, self.temp_isometric_state, self.temp_gamma))
                ion_index = [i for i, n in enumerate(self.spectrum_source.data['ion']) if n == self.temp_ion]
                iso_index = [i for i, n in enumerate(self.spectrum_source.data['isometric']) if n == self.temp_isometric_state]
                index = np.intersect1d(ion_index, iso_index)
                self.ion_harmonics.data = {'xs': [self.spectrum_source.data['xs'][_index] for _index in index], 'ys': [self.spectrum_source.data['ys'][_index] for _index in index], 'ion': [self.spectrum_source.data['ion'][_index] for _index in index], 'isometric': [self.spectrum_source.data['isometric'][_index] for _index in index], 'peak_loc': [self.spectrum_source.data['peak_loc'][_index] for _index in index], 'weight': [self.spectrum_source.data['weight'][_index] for _index in index], 'yield': [self.spectrum_source.data['yield'][_index] for _index in index], 'harmonic': [self.spectrum_source.data['harmonic'][_index] for _index in index], 'rev_freq': [self.spectrum_source.data['rev_freq'][_index] for _index in index], 'half_life': [self.spectrum_source.data['half_life'][_index] for _index in index]}
            except:
                pass

        self.spectrum_source = ColumnDataSource(data=self._wrap_data(1))
        self.ion_harmonics = ColumnDataSource(data=self._wrap_data(-1))
        ion_tooltip = [
                ("ion", '@ion'+'('+'@isometric'+')'),
                ("peak location", '@peak_loc'+' kHz'),
                ("weight", '@weight'),
                ("yield", '@yield'),
                ("harmonic", '@harmonic'),
                ("revolution frequency", '@rev_freq' + ' MHz'),
                ("half life", '@half_life')
        ]
        self.p_spectrum_default_log = figure(width=1000, height=300, title='Simulated Spectrum (lifetime > 10 ms)', tools='pan, crosshair, tap, box_zoom, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save, hover', x_range=(-self.iid.span/2,self.iid.span/2), y_axis_type='log', output_backend='webgl')
        self.p_spectrum_default_log.tools[-1].tooltips=ion_tooltip            
        self.p_spectrum_default_log.tools[-1].attachment = 'vertical'
        self.p_spectrum_default_log.tools[-1].point_policy = 'follow_mouse'
        self.p_spectrum_default_log.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)
        self.p_spectrum_default_log.yaxis.axis_label = "psd [arb. unit]"
        self.p_spectrum_default_log.patches(xs='xs', ys='ys', hover_color='darkorange', selection_color='red', source=self.spectrum_source, color="dimgray")
        self.p_spectrum_default_log.patches(xs='xs', ys='ys', source=self.ion_harmonics, color='goldenrod')

        self.p_spectrum_default_linear = figure(width=1000, height=300, title='Simulated Spectrum (lifetime > 10 ms)', tools='pan, crosshair, tap, box_zoom, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save, hover', x_range=(-self.iid.span/2,self.iid.span/2), output_backend='webgl')
        self.p_spectrum_default_linear.tools[-1].tooltips=ion_tooltip            
        self.p_spectrum_default_linear.tools[-1].attachment = 'vertical'
        self.p_spectrum_default_linear.tools[-1].point_policy = 'follow_mouse'
        self.p_spectrum_default_linear.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)
        self.p_spectrum_default_linear.yaxis.axis_label = "psd [arb. unit]"
        self.p_spectrum_default_linear.patches(xs='xs', ys='ys', hover_color='darkorange', selection_color='red', source=self.spectrum_source, color="dimgray")
        self.p_spectrum_default_linear.patches(xs='xs', ys='ys', source=self.ion_harmonics, color='goldenrod')

        self.labels_default = LabelSet(x='peak_loc', y='peak_max', source=self.spectrum_source, text='ion_label', x_offset=0, y_offset=0)
        self.p_spectrum_default_log.add_layout(self.labels_default)
        self.p_spectrum_default_linear.add_layout(self.labels_default)

        self.p_yield_default = figure(width=550, height=550, title='Ion Yield', tools='pan, box_zoom, tap, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save', x_range=(-0.5,177.5), y_range=(-0.5,118.5), aspect_ratio=1., tooltips=ion_tooltip, output_backend='webgl')
        self.p_yield_default.rect(x='N', y='Z', fill_color='color', source=self.spectrum_source, line_color='lightgray', width=1., height=1.)
        try:
            yield_top = int(np.log10(np.max(self.spectrum_source.data['total_yield'])))
        except:
            yield_top = 1
        self.default_colorBar = LogColorMapper(palette=YlOrRd9, high=10**(yield_top-8), low=10**(yield_top+1))
        color_bar_default = ColorBar(color_mapper=self.default_colorBar)
        self.p_yield_default.add_layout(color_bar_default, "left")

        self.spectrum_source.selected.on_change("indices", selected_ion)

        columns = [
                TableColumn(field='ion', title='ion'),
                TableColumn(field='isometric', title='isometric state'),
                TableColumn(field='half_life', title='half life'),
                TableColumn(field='weight', title='weight'),
                TableColumn(field='yield', title='yield'),
                TableColumn(field='harmonic', title='harmonic'),
                TableColumn(field='peak_loc', title='peak loc [kHz]'),
                TableColumn(field='rev_freq', title='rev freq [MHz]')
        ]
        self.p_table_default = DataTable(source=self.spectrum_source, columns=columns, width=1000, height=300, frozen_columns=3, index_position=-1, sortable=True, selectable=True, stylesheets=[InlineStyleSheet(css='.slick-cell.selected {background-color: #F1B6B9;}')])
        # inital tabpanel
        self.tabpanel_ec_off = TabPanel(child=column([self.p_spectrum_default_linear, self.p_spectrum_default_log, self.p_table_default]), title='EC OFF')
        self.tabpanel_yield_ec_off = TabPanel(child=self.p_yield_default, title='EC OFF')

    def _initial_RevTime_spectrum(self):
        self.TOF_source = ColumnDataSource(data=self._wrap_data(-2))
        ion_tooltip = [
                ("ion", '@ion'+'('+'@isometric'+')'),
                ("yield", '@yield'),
                ("revolution time", '@rev_time' + ' ns'),
                ("half life", '@half_life')
        ]
        self.p_spectrum_TOF_log = figure(width=1000, height=300, title='Simulated Spectrum (lifetime > 10 ms)', tools='pan, crosshair, tap, box_zoom, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save, hover', y_axis_type='log', output_backend='webgl')
        self.p_spectrum_TOF_log.tools[-1].tooltips=ion_tooltip            
        self.p_spectrum_TOF_log.tools[-1].attachment = 'vertical'
        self.p_spectrum_TOF_log.tools[-1].point_policy = 'follow_mouse'
        self.p_spectrum_TOF_log.xaxis.axis_label = "revolution time [ns]"
        self.p_spectrum_TOF_log.yaxis.axis_label = "counts/s"
        self.p_spectrum_TOF_log.patches(xs='xs', ys='ys', hover_color='darkorange', selection_color='red', source=self.TOF_source, color="dimgray")

        self.p_spectrum_TOF_linear = figure(width=1000, height=300, title='Simulated Spectrum (lifetime > 10 ms)', tools='pan, crosshair, tap, box_zoom, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save, hover', output_backend='webgl')
        self.p_spectrum_TOF_linear.tools[-1].tooltips=ion_tooltip            
        self.p_spectrum_TOF_linear.tools[-1].attachment = 'vertical'
        self.p_spectrum_TOF_linear.tools[-1].point_policy = 'follow_mouse'
        self.p_spectrum_TOF_linear.xaxis.axis_label = "revolution time [ns]"
        self.p_spectrum_TOF_linear.yaxis.axis_label = "counts/s"
        self.p_spectrum_TOF_linear.patches(xs='xs', ys='ys', hover_color='darkorange', selection_color='red', source=self.TOF_source, color="dimgray")

        self.labels_TOF = LabelSet(x='rev_time', y='peak_max', source=self.TOF_source, text='ion_label', x_offset=0, y_offset=0)
        self.p_spectrum_TOF_log.add_layout(self.labels_TOF)
        self.p_spectrum_TOF_linear.add_layout(self.labels_TOF)

        self.p_yield_TOF = figure(width=550, height=550, title='Ion Yield', tools='pan, box_zoom, tap, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save', x_range=(-0.5,177.5), y_range=(-0.5,118.5), aspect_ratio=1., tooltips=ion_tooltip, output_backend='webgl')
        self.p_yield_TOF.rect(x='N', y='Z', fill_color='color', source=self.TOF_source, line_color='lightgray', width=1., height=1.)
        try:
            yield_top = int(np.log10(np.max(self.TOF_source.data['total_yield'])))
        except:
            yield_top = 1
        self.TOF_colorBar = LogColorMapper(palette=YlOrRd9, high=10**(yield_top-8), low=10**(yield_top+1))
        color_bar_TOF = ColorBar(color_mapper=self.TOF_colorBar)
        self.p_yield_TOF.add_layout(color_bar_TOF, "left")

        #self.spectrum_source.selected.on_change("indices", selected_ion)

        columns = [
                TableColumn(field='ion', title='ion'),
                TableColumn(field='isometric', title='isometric state'),
                TableColumn(field='half_life', title='half life'),
                TableColumn(field='yield', title='yield'),
                TableColumn(field='rev_time', title='rev time [ns]')
        ]
        self.p_table_TOF = DataTable(source=self.TOF_source, columns=columns, width=1000, height=300, frozen_columns=3, index_position=-1, sortable=True, selectable=True, stylesheets=[InlineStyleSheet(css='.slick-cell.selected {background-color: #F1B6B9;}')])
        # inital tabpanel
        self.tabpanel_TOF = TabPanel(child=column([self.p_spectrum_TOF_linear, self.p_spectrum_TOF_log, self.p_table_TOF]), title='TOF')
        self.tabpanel_yield_TOF = TabPanel(child=self.p_yield_TOF, title='TOF')

    def _update_spectrum_labels(self):
        self.p_spectrum_default_log.x_range.start = -self.iid.span/2
        self.p_spectrum_default_log.x_range.end = self.iid.span/2
        self.p_spectrum_default_log.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)
        self.p_spectrum_default_linear.x_range.start = -self.iid.span/2
        self.p_spectrum_default_linear.x_range.end = self.iid.span/2
        self.p_spectrum_default_linear.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)
        self.p_spectrum_cooler_log.x_range.start = -self.iid.span/2
        self.p_spectrum_cooler_log.x_range.end = self.iid.span/2
        self.p_spectrum_cooler_log.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)
        self.p_spectrum_cooler_linear.x_range.start = -self.iid.span/2
        self.p_spectrum_cooler_linear.x_range.end = self.iid.span/2
        self.p_spectrum_cooler_linear.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)

    def _initial_ec_on(self):
        def selected_ion(attr, old, new):
            try:
                self.temp_ion = self.cooler_source.data['ion'][new[0]]
                self.temp_isometric_state = self.cooler_source.data['isometric'][new[0]]
                self.temp_harmonic = self.cooler_source.data['harmonic'][new[0]]
                self.temp_gamma = self.cooler_source.data['pseudo_gamma'][new[0]]
                self.input_gamma_setting.value = self.temp_gamma
                print("{:}({:}), γ: {:.3f}".format(self.temp_ion, self.temp_isometric_state, self.temp_gamma))
                self._log("Ion Selected: {:}({:}), γ: {:.3f}".format(self.temp_ion, self.temp_isometric_state, self.temp_gamma))
                ion_index = [i for i, n in enumerate(self.cooler_source.data['ion']) if n == self.temp_ion]
                iso_index = [i for i, n in enumerate(self.cooler_source.data['isometric']) if n == self.temp_isometric_state]
                index = np.intersect1d(ion_index, iso_index)
                self.cooler_harmonics.data = {'xs': [self.cooler_source.data['xs'][_index] for _index in index], 'ys': [self.cooler_source.data['ys'][_index] for _index in index], 'ion': [self.cooler_source.data['ion'][_index] for _index in index], 'isometric': [self.cooler_source.data['isometric'][_index] for _index in index], 'peak_loc': [self.cooler_source.data['peak_loc'][_index] for _index in index], 'weight': [self.cooler_source.data['weight'][_index] for _index in index], 'yield': [self.cooler_source.data['yield'][_index] for _index in index], 'harmonic': [self.cooler_source.data['harmonic'][_index] for _index in index], 'rev_freq': [self.cooler_source.data['rev_freq'][_index] for _index in index], 'half_life': [self.cooler_source.data['half_life'][_index] for _index in index]}
            except:
                pass

        self.cooler_source = ColumnDataSource(data=self._wrap_data(0))
        self.cooler_harmonics = ColumnDataSource(data=self._wrap_data(-1))
        ion_tooltip = [
                ("ion", '@ion'+'('+'@isometric'+')'),
                ("peak location", '@peak_loc'+' kHz'),
                ("weight", '@weight'),
                ("yield", '@yield'),
                ("harmonic", '@harmonic'),
                ("revolution frequency", '@rev_freq' + ' MHz'),
                ("half life", '@half_life')
        ]
        self.p_spectrum_cooler_log = figure(width=1000, height=300, title='Simulated Spectrum (lifetime > 1 sec)', tools='pan, tap, box_zoom, crosshair, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save, hover', x_range=(-self.iid.span/2,self.iid.span/2), y_axis_type='log', output_backend='webgl')
        self.p_spectrum_cooler_log.tools[-1].tooltips=ion_tooltip            
        self.p_spectrum_cooler_log.tools[-1].attachment = 'vertical'
        self.p_spectrum_cooler_log.tools[-1].point_policy = 'follow_mouse'
        self.p_spectrum_cooler_log.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)
        self.p_spectrum_cooler_log.yaxis.axis_label = "psd [arb. unit]"
        self.p_spectrum_cooler_log.patches(xs='xs', ys='ys', hover_color='darkorange', selection_color='lime', source=self.cooler_source, color="deepskyblue")
        self.p_spectrum_cooler_log.patches(xs='xs', ys='ys', source=self.cooler_harmonics, color='goldenrod')
        self.p_spectrum_cooler_linear = figure(width=1000, height=300, title='Simulated Spectrum (lifetime > 1 sec)', tools='pan, tap, crosshair, box_zoom, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save, hover', x_range=(-self.iid.span/2,self.iid.span/2), output_backend='webgl')
        self.p_spectrum_cooler_linear.tools[-1].tooltips=ion_tooltip            
        self.p_spectrum_cooler_linear.tools[-1].attachment = 'vertical'
        self.p_spectrum_cooler_linear.tools[-1].point_policy = 'follow_mouse'
        self.p_spectrum_cooler_linear.xaxis.axis_label = "{:} MHz [kHz]".format(self.iid.cen_freq)
        self.p_spectrum_cooler_linear.yaxis.axis_label = "psd [arb. unit]"
        self.p_spectrum_cooler_linear.patches(xs='xs', ys='ys', hover_color='darkorange', selection_color='lime', source=self.cooler_source, color="deepskyblue")
        self.p_spectrum_cooler_linear.patches(xs='xs', ys='ys', source=self.cooler_harmonics, color='goldenrod')

        self.labels_cooler = LabelSet(x='peak_loc', y='peak_max', source=self.cooler_source, text='ion_label', x_offset=0, y_offset=0)
        self.p_spectrum_cooler_log.add_layout(self.labels_cooler)
        self.p_spectrum_cooler_linear.add_layout(self.labels_cooler)

        self.p_yield_cooler = figure(width=550, height=550, title='Ion Yield', tools='pan, box_zoom, tap, wheel_zoom, zoom_in, zoom_out, undo, redo, reset, save', x_range=(-0.5,177.5), y_range=(-0.5,118.5), aspect_ratio=1., tooltips=ion_tooltip, output_backend='webgl')
        self.p_yield_cooler.rect(x='N', y='Z', fill_color='color', source=self.cooler_source, line_color='lightgray', width=1., height=1.)
        try:
            yield_top = int(np.log10(np.max(self.cooler_source.data['total_yield'])))
        except:
            yield_top = 1
        self.cooler_colorBar = LogColorMapper(palette=YlGnBu9, high=10**(yield_top-8), low=10**(yield_top+1))
        color_bar_cooler = ColorBar(color_mapper=self.cooler_colorBar)
        self.p_yield_cooler.add_layout(color_bar_cooler, "left")
        
        self.cooler_source.selected.on_change("indices", selected_ion)

        columns = [
                TableColumn(field='ion', title='ion'),
                TableColumn(field='isometric', title='isometric state'),
                TableColumn(field='half_life', title='half life'),
                TableColumn(field='weight', title='weight'),
                TableColumn(field='yield', title='yield'),
                TableColumn(field='harmonic', title='harmonic'),
                TableColumn(field='peak_loc', title='peak loc [kHz]'),
                TableColumn(field='rev_freq', title='rev freq [MHz]')
        ]
        self.cooler_source.selected.on_change("indices", selected_ion)
        self.p_table_cooler = DataTable(source=self.cooler_source, columns=columns, width=1000, height=300, frozen_columns=3, index_position=-1, sortable=True, selectable=True, stylesheets=[InlineStyleSheet(css='.slick-row.odd {background-color: #C5E6F9;} .slick-cell.selected {background-color: #CFF9A8;}')])
        # inital tabpanel
        self.tabpanel_ec_on = TabPanel(child=column([self.p_spectrum_cooler_linear, self.p_spectrum_cooler_log, self.p_table_cooler]), title='EC ON')
        self.tabpanel_yield_ec_on = TabPanel(child=self.p_yield_cooler, title='EC ON')

    def _update(self, data_type=1):
        if data_type == 1:
            data = self._wrap_data(1)
            self.spectrum_source.data = data
            self.p_table_default.source.data = data
            yield_top = int(np.log10(np.max(data['total_yield'])))
            self.default_colorBar.low = 10**(yield_top+1)
            self.default_colorBar.high = 10**(yield_top-8)
            self.ion_harmonics.data = self._wrap_data(-1)
        elif data_type == 0:
            data = self._wrap_data(0)
            self.cooler_source.data = data 
            self.p_table_cooler.source.data = data
            yield_top = int(np.log10(np.max(data['total_yield'])))
            self.cooler_colorBar.low = 10**(yield_top+1)
            self.cooler_colorBar.high = 10**(yield_top-8)
            self.cooler_harmonics.data = self._wrap_data(-1)
        else:
            data = self._wrap_data(-2)
            self.TOF_source.data = data
            self.p_table_TOF.source.data = data
            yield_top = int(np.log10(np.max(data['total_yield'])))
            self.TOF_colorBar.low = 10**(yield_top+1)
            self.TOF_colorBar.high = 10**(yield_top-8)
    
    def panel_control(self):
        '''
        panel to control the bokeh show
        control panel for file
        Brho, peak location calibrate
        '''
        # input for Brho
        self.checkbox_Brho_input = Checkbox(label='Using Bρ for calibrate', height=20, active=False)
        self.input_Brho = NumericInput(value=self.iid.Brho, height=50, low=1., high=15., mode='float', title='Bρ [Tm]', disabled=True)
        def update_Brho_log(attr, old, new):
            self._log('calibrate ...')
        def update_Brho(attr, old, new):
            print('calibrate Brho...')
            self.iid.calibrate_Brho(float(new), self.checkbox_ec_on.active)
            self._update(1)
            if self.checkbox_ec_on.active:
                self._update(0)
            print('calibrate complete!')
            self._log('calibrate complete!')
        self.input_Brho.on_change('value', update_Brho_log, update_Brho)
        def set_Brho_input(attr, old, new):
            if self.checkbox_Brho_input.active:
                self.input_Brho.disabled = False
                self.input_peakloc.disabled = True
                self.button_calibrate.disabled = True
            else:
                self.input_Brho.disabled = True
                self.input_peakloc.disabled = False
                self.button_calibrate.disabled = False
        self.checkbox_Brho_input.on_change('active', set_Brho_input)

        # button for calibrate of ion
        self.input_peakloc = NumericInput(value=0, height=50, low=-self.iid.span/2, high=self.iid.span, mode='float', title='peak location [kHz]')
        self.button_calibrate = Button(label='calibrate', height=50, width=80, button_type='primary')
        def calibrate_ion():
            try:
                print('calibrate ion peak loc...')
                Brho = self.iid.calibrate_peak_loc(self.temp_ion, self.temp_isometric_state, self.input_peakloc.value, self.temp_harmonic)
                self.input_Brho.value = Brho
            except:
                print('warning: no ion selected')
                self._log('warning: select one ion for peak location calibrate first!')
        self.button_calibrate.on_event(ButtonClick, calibrate_ion)

        # button for ecooler
        self.checkbox_ec_on = Checkbox(label='EC on', height=20, active=False)
        self.input_gamma_setting = NumericInput(value=self.iid.gamma_setting, height=50, low=1.0001, high=5.0000, mode='float', title='γ setting', disabled=True)
        m_over_q = self.iid.Brho / self.iid.gamma_setting / np.sqrt(1 - 1/self.iid.gamma_setting**2) / self.iid.c / self.iid.u2kg * self.iid.e
        self.input_m_over_q = NumericInput(value=m_over_q, height=50, low=0., high=100., mode='float', title='m/q', disabled=True)
        self.input_delta_v_over_v = NumericInput(value=self.iid.delta_v_over_v, height=50, low=1e-8, high=1e-5, mode='float', title='Δv/v', disabled=True)
        self.button_set = Button(label='set', height=50, width=80, button_type='primary', disabled=True)
        def set_velocity():
            if self.checkbox_ec_on.active:
                print('calibrate Δv/v ...')
                self.iid.calibrate_ecooler(self.input_gamma_setting.value, self.input_delta_v_over_v.value)
                self._update(0)
                print('calibrate complete!')
                self._log('setting complete!')
        self.button_set.on_event(ButtonClick, set_velocity)
        def set_gamma_setting(attr, old, new):
            self.input_m_over_q.value = self.input_Brho.value / float(new) / np.sqrt(1 - 1/float(new)**2) / self.iid.c / self.iid.u2kg * self.iid.e
        self.input_gamma_setting.on_change('value', set_gamma_setting)
        def set_m_over_q(attr, old, new):
            self.input_gamma_setting.value = np.sqrt(1 + (self.input_Brho.value / float(new) / self.iid.c / self.iid.u2kg * self.iid.e)**2)
        self.input_m_over_q.on_change('value', set_m_over_q)
        def set_ec_on(attr, old, new):
            if self.checkbox_ec_on.active:
                self.input_gamma_setting.disabled = False
                self.input_m_over_q.disabled = False
                self.input_delta_v_over_v.disabled = False
                self.button_set.disabled = False
            else:
                self.input_gamma_setting.disabled = True
                self.input_m_over_q.disabled = True
                self.input_delta_v_over_v.disabled = True
                self.button_set.disabled = True
        self.checkbox_ec_on.on_change('active', set_ec_on)

        # button for global setting
        self.input_gamma_t = NumericInput(value=self.iid.gamma_t, height=50, low=1.0001, high=5.0000, mode='float', title='γt')
        self.input_alpha_p = NumericInput(value=1/self.iid.gamma_t**2, height=50, low=0.04, high=0.99999, mode='float', title='αp')
        self.input_delta_Brho_over_Brho = NumericInput(value=self.iid.delta_Brho_over_Brho, height=50, low=0.01, high=10.00, mode='float', title='ΔΒρ/Βρ, %')
        self.checkbox_log_or_linear = Checkbox(label='log scale', height=20, active=True)
        def update_gamma_t(attr, old, new):
            self.input_alpha_p.value = 1/float(new)**2
            print('update γt ...')
            self.iid.update_gamma_t(float(new), self.checkbox_ec_on.active)
            self._update(1)
            self._update(-2)
            if self.checkbox_ec_on.active:
                self._update(0)
            print('update complete!')
            self._log('update complete!')
        self.input_gamma_t.on_change('value', update_gamma_t)
        def update_alpha_p(attr, old, new):
            self.input_gamma_t.value = 1/np.sqrt(float(new))
        self.input_alpha_p.on_change('value', update_alpha_p)
        def update_delta_Brho_over_Brho(attr, old, new):
            print('update ΔΒρ/Βρ ...')
            self.iid.update_delta_Brho_over_Brho(float(new), self.checkbox_ec_on.active)
            self._update(1)
            self._update(-2)
            if self.checkbox_ec_on.active:
                self._update(0)
            print('update complete!')
            self._log('update complete!')
        self.input_delta_Brho_over_Brho.on_change('value', update_delta_Brho_over_Brho)
        def set_log(attr, old, new):
            if self.checkbox_log_or_linear.active:
                self.p_spectrum_default_log.visible = True
                self.p_spectrum_default_linear.visible = False
                self.p_spectrum_cooler_log.visible = True
                self.p_spectrum_cooler_linear.visible = False
                self.p_spectrum_TOF_log.visible = True
                self.p_spectrum_TOF_linear.visible = False
            else:
                self.p_spectrum_default_log.visible = False
                self.p_spectrum_default_linear.visible = True
                self.p_spectrum_cooler_log.visible = False
                self.p_spectrum_cooler_linear.visible = True
                self.p_spectrum_TOF_log.visible = False
                self.p_spectrum_TOF_linear.visible = True
        self.checkbox_log_or_linear.on_change('active', set_log)

        self.input_cen_freq = NumericInput(value=self.iid.cen_freq, height=50, low=20, high=450, mode='float', title='center frequency [MHz]')
        self.input_span = NumericInput(value=self.iid.span, height=50, low=10, high=20000, mode='float', title='span [kHz]')
        self.input_L_CSRe = NumericInput(value=self.iid.L_CSRe, height=50, low=10, high=400, mode='float', title='length of Ring [m]')
        def update_cen_freq(attr, old, new):
            print('update center frequency ...')
            self.iid.update_cen_freq(float(new), self.checkbox_ec_on.active)
            self._update_spectrum_labels()
            self._update(1)
            if self.checkbox_ec_on.active:
                self._update(0)
            print('update complete!')
            self._log('update complete!')
        self.input_cen_freq.on_change('value', update_cen_freq)
        def update_span(attr, old, new):
            print('update span ...')
            self.iid.update_span(float(new), self.checkbox_ec_on.active)
            self._update_spectrum_labels()
            self._update(1)
            if self.checkbox_ec_on.active:
                self._update(0)
            print('update complete!')
            self._log('update complete!')
        self.input_span.on_change('value', update_span)
        def update_L_CSRe(attr, old, new):
            print('update length of Ring ...')
            self.iid.update_L_CSRe(float(new), self.checkbox_ec_on.active)
            self._update_spectrum_labels()
            self._update(1)
            self._update(-2)
            if self.checkbox_ec_on.active:
                self._update(0)
            print('update complete!')
        self.input_L_CSRe.on_change('value', update_L_CSRe)

        result = self.iid.cur.execute("SELECT DISTINCT ION, ISOMERIC FROM OBSERVEDION").fetchall()
        ion_completion = ["{:}({:})".format(ion, isometric) for ion, isometric in result]
        self.input_ion = AutocompleteInput(completions=ion_completion, title='ion')
        self.checkbox_TOF_on = Checkbox(label='TOF on', height=25, active=False)
        self.button_find_ion = Button(label='find', height=50, width=80, button_type='primary')
        self.button_reset_ion = Button(label='reset', height=50, width=80, button_type='primary')
        def find_ion():
            if self.input_ion.value != '':
                ion, isometric = self.input_ion.value.split('(')
                print('{:}({:})'.format(ion, isometric[:-1]))
                if self.checkbox_TOF_on:
                    ion_index = [i for i, n in enumerate(self.TOF_source.data['ion']) if n == ion]
                    iso_index = [i for i, n in enumerate(self.TOF_source.data['isometric']) if n == isometric[:-1]]
                    index = np.intersect1d(ion_index, iso_index)
                    if len(index) < 1:
                        self._log('no ion available in the TOF spectrum!')
                        return
                    self.TOF_source.selected.indices = [index[0]]
                    return
                if self.checkbox_ec_on.active:
                    ion_index = [i for i, n in enumerate(self.cooler_source.data['ion']) if n == ion]
                    iso_index = [i for i, n in enumerate(self.cooler_source.data['isometric']) if n == isometric[:-1]]
                    index = np.intersect1d(ion_index, iso_index)
                    if len(index) < 1:
                        self._log('no ion available in the Schottky spectrum (EC on)!')
                        return
                    self.cooler_harmonics.data = {'xs': [self.cooler_source.data['xs'][_index] for _index in index], 'ys': [self.cooler_source.data['ys'][_index] for _index in index], 'ion': [self.cooler_source.data['ion'][_index] for _index in index], 'isometric': [self.cooler_source.data['isometric'][_index] for _index in index], 'peak_loc': [self.cooler_source.data['peak_loc'][_index] for _index in index], 'weight': [self.cooler_source.data['weight'][_index] for _index in index], 'yield': [self.cooler_source.data['yield'][_index] for _index in index], 'harmonic': [self.cooler_source.data['harmonic'][_index] for _index in index], 'rev_freq': [self.cooler_source.data['rev_freq'][_index] for _index in index], 'half_life': [self.cooler_source.data['half_life'][_index] for _index in index]}
                    self.cooler_source.selected.indices = [index[0]]
                else:
                    ion_index = [i for i, n in enumerate(self.spectrum_source.data['ion']) if n == ion]
                    iso_index = [i for i, n in enumerate(self.spectrum_source.data['isometric']) if n == isometric[:-1]]
                    index = np.intersect1d(ion_index, iso_index)
                    if len(index) < 1:
                        self._log('no ion available in the Schottky spectrum (EC off)!')
                        return
                    self.ion_harmonics.data = {'xs': [self.spectrum_source.data['xs'][_index] for _index in index], 'ys': [self.spectrum_source.data['ys'][_index] for _index in index], 'ion': [self.spectrum_source.data['ion'][_index] for _index in index], 'isometric': [self.spectrum_source.data['isometric'][_index] for _index in index], 'peak_loc': [self.spectrum_source.data['peak_loc'][_index] for _index in index], 'weight': [self.spectrum_source.data['weight'][_index] for _index in index], 'yield': [self.spectrum_source.data['yield'][_index] for _index in index], 'harmonic': [self.spectrum_source.data['harmonic'][_index] for _index in index], 'rev_freq': [self.spectrum_source.data['rev_freq'][_index] for _index in index], 'half_life': [self.spectrum_source.data['half_life'][_index] for _index in index]}
                    self.spectrum_source.selected.indices = [index[0]]
        self.button_find_ion.on_event(ButtonClick, find_ion)
        def reset_ion():
            self.ion_harmonics.data = self._wrap_data(-1)
            self.cooler_harmonics.data = self._wrap_data(-1)
        self.button_reset_ion.on_event(ButtonClick, reset_ion)
        def set_TOF_on(attr, old, new):
            if self.checkbox_TOF_on.active:
                self.checkbox_Brho_input.active = True
                self.input_cen_freq.disabled = True
                self.input_span.disabled = True
            else:
                self.checkbox_Brho_input.active = False
                self.input_cen_freq.disabled = False
                self.input_span.disabled = False
        self.checkbox_TOF_on.on_change('active', set_TOF_on)

        #result = self.iid.cur.execute("SELECT min(PEAKMAX) FROM ISOCHRONOUSION").fetchone()[0]
        self.input_show_threshold = NumericInput(value=1e-16, low=1e-16, high=1e16, height=50, mode='float', title='threshold')
        self.checkbox_labels_on = Checkbox(label='show labels', height=25, active=False)
        def set_threshold(attr, old, new):
            self._update(1)
            self._update(-2)
            if self.checkbox_ec_on.active:
                self._update(0)
        self.input_show_threshold.on_change('value', set_threshold)
        def set_labels_on(attr, old, new):
            if self.checkbox_labels_on.active:
                self.labels_default.visible = True
                self.labels_cooler.visible = True
                self.labels_TOF.visible = True
            else:
                self.labels_default.visible = False
                self.labels_cooler.visible = False
                self.labels_TOF.visible = False
        self.checkbox_labels_on.on_change('active', set_labels_on)
        
        self.div_log = Div(text='', width=300, height=50, background='darkorange')

if __name__ == '__main__':
    curdoc().add_root(Bokeh_show('./Test_CSRe_173Er67.lpp', 243., 3000, 1.34, 0.2, 1.34, 1e-6)._show())
